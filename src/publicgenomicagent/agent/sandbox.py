"""Sandbox consumer: read synthetic case catalogs.

Sandbox is a separate project that generates synthetic genomic
cases. PublicGenomicAgent consumes them as external inputs.

This module implements the consumer side of the contract:

  - SandboxCatalog: lists and indexes cases.
  - SandboxCase: represents a single case.
  - build_case_manifest: converts a SandboxCase to a CaseManifest
    ready to feed the AgentLoop.

The generator is opaque: we do not know (or need to know) how
the BAMs and manifests were produced. We only read them.

Contract reference: docs/sandbox_contract.md (or the Sandbox
repository's SANDBOX_CONTRACT.md).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from ..tools.base import CaseManifest, GenomicRange, IndividualSpec


# Zygosity codes used by Sandbox -> genotype strings in VCF
ZYGOSITY_TO_GT = {
    "HOM_WT": "0/0",
    "HET": "0/1",
    "HOM_ALT": "1/1",
    "HOM_DEL": "1/1",
}


class SandboxError(RuntimeError):
    """Error reading or validating a Sandbox case."""


@dataclass
class SandboxCase:
    """A single case loaded from a Sandbox catalog.

    Attributes:
      case_id:    deterministic ID = f"{gene}:{spdi}"
      gene:       HGNC symbol (e.g. NPHP1)
      spdi:       canonical SPDI of the key variant
      title:      human-readable title
      clinical_text: narrative clinical description
      chrom:      chromosome used in the synthetic reference
      zygosity:   dict sample -> Sandbox zygosity code
      hpo_expected: list of {code, term}
      path:       absolute path to the case directory
      bams:       dict sample -> absolute BAM path
    """
    case_id: str
    gene: str
    spdi: str
    title: str
    clinical_text: str
    chrom: str
    zygosity: dict[str, str]
    hpo_expected: list[dict[str, str]]
    path: Path
    bams: dict[str, Path] = field(default_factory=dict)

    def expected_genotypes(self) -> dict[str, str]:
        """Sandbox zygosity -> VCF GT strings."""
        return {
            sample: ZYGOSITY_TO_GT.get(zyg, "?")
            for sample, zyg in self.zygosity.items()
        }

    def expected_proband_gt(self) -> str | None:
        """Expected proband genotype (VCF string)."""
        if not self.zygosity:
            return None
        # Proband is conventionally the last sample or "proband".
        if "proband" in self.zygosity:
            return ZYGOSITY_TO_GT.get(self.zygosity["proband"])
        # Fallback: pick the affected genotype
        for sample, zyg in self.zygosity.items():
            if zyg in ("HOM_DEL", "HOM_ALT"):
                return ZYGOSITY_TO_GT.get(zyg)
        return None

    def sandbox_zygosity(self, sample: str) -> str | None:
        return self.zygosity.get(sample)


def _parse_spdi(spdi: str) -> tuple[str, int, int, str]:
    """Parse SPDI: RefSeq:pos:deleted:inserted.

    Returns (refseq, pos, del_len, ins_seq).
    dele_len = len(deleted) if numeric, else len(string).
    """
    parts = spdi.split(":")
    if len(parts) < 4:
        raise SandboxError(f"Invalid SPDI: {spdi}")
    refseq = parts[0]
    try:
        pos = int(parts[1])
    except ValueError as e:
        raise SandboxError(f"Invalid SPDI position: {parts[1]}") from e
    del_field = parts[2]
    if del_field.isdigit():
        del_len = int(del_field)
    else:
        del_len = len(del_field)
    ins_seq = parts[3]
    return refseq, pos, del_len, ins_seq


def _roi_from_ground_truth(ground_truth: dict, *, padding: int = 100) -> GenomicRange | None:
    """Build a GenomicRange from the Sandbox ground_truth.

    Uses chrom + SPDI position, padded by `padding` bp on each side.
    For large deletions the range covers the whole deletion.
    """
    chrom = ground_truth.get("chrom")
    if not chrom:
        return None
    spdi = ground_truth.get("spdi") or ground_truth.get("variant_spdi")
    # SPDI may be at top level of the manifest, not inside ground_truth
    return None  # filled by caller with top-level spdi


class SandboxCatalog:
    """Read a Sandbox catalog from a root directory.

    Prefers `registry.json` (fast index built by Sandbox).
    Falls back to walking `genes/*/*/manifest.yaml` if the
    registry is missing or stale.
    """

    def __init__(self, root: Path | str):
        self.root = Path(root).expanduser().resolve()
        if not self.root.exists():
            raise SandboxError(f"Sandbox root does not exist: {self.root}")

    def list_cases(self) -> list[SandboxCase]:
        """Return all cases in the catalog, sorted by gene + spdi."""
        registry = self.root / "registry.json"
        if registry.exists():
            cases = self._load_from_registry(registry)
        else:
            cases = self._load_from_filesystem()
        cases.sort(key=lambda c: (c.gene, c.spdi))
        return cases

    def get_case(self, case_id: str) -> SandboxCase:
        """Look up a case by case_id ("gene:spdi") or gene name.

        If multiple cases match a gene name, raises. Use the
        full case_id ("gene:spdi") to disambiguate.
        """
        cases = self.list_cases()
        # Exact case_id match
        for c in cases:
            if c.case_id == case_id:
                return c
        # Gene-name match (must be unique)
        matches = [c for c in cases if c.gene == case_id]
        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise SandboxError(f"No case found for: {case_id}")
        raise SandboxError(
            f"Multiple cases for gene {case_id}. Use full ID: "
            + ", ".join(c.case_id for c in matches)
        )

    # --- internals --------------------------------------------------------

    def _load_from_registry(self, registry_path: Path) -> list[SandboxCase]:
        data = json.loads(registry_path.read_text(encoding="utf-8"))
        cases: list[SandboxCase] = []
        for gene, entries in data.items():
            for entry in entries:
                try:
                    case = self._case_from_entry(gene, entry)
                except Exception as e:  # noqa: BLE001
                    raise SandboxError(
                        f"Error parsing entry for gene {gene}: {e}"
                    ) from e
                cases.append(case)
        return cases

    def _case_from_entry(self, gene: str, entry: dict) -> SandboxCase:
        spdi = entry.get("spdi", "")
        if not spdi:
            raise SandboxError(f"Entry without spdi for gene {gene}")
        case_id = f"{gene}:{spdi}"

        # Paths
        paths = entry.get("paths", {})
        folder = paths.get("folder")
        if not folder:
            raise SandboxError(f"Entry without paths.folder: {case_id}")
        case_path = (self.root / folder).resolve()
        if not case_path.exists():
            raise SandboxError(f"Case folder missing: {case_path}")

        bams: dict[str, Path] = {}
        for sample in ("father", "mother", "proband"):
            key = f"{sample}_bam"
            rel = paths.get(key)
            if rel:
                p = (self.root / rel).resolve()
                if p.exists():
                    bams[sample] = p

        # Ground truth
        gt = entry.get("ground_truth", {})
        chrom = gt.get("chrom", "")
        zygosity = gt.get("zygosity", {})

        # HPO
        hpo = entry.get("hpo_expected", []) or []

        return SandboxCase(
            case_id=case_id,
            gene=gene,
            spdi=spdi,
            title=entry.get("title", ""),
            clinical_text=entry.get("clinical_text", ""),
            chrom=chrom,
            zygosity=zygosity,
            hpo_expected=hpo,
            path=case_path,
            bams=bams,
        )

    def _load_from_filesystem(self) -> list[SandboxCase]:
        """Fallback: walk genes/*/*/manifest.yaml."""
        genes_dir = self.root / "genes"
        if not genes_dir.exists():
            raise SandboxError(
                f"No registry.json and no genes/ folder in {self.root}"
            )
        cases: list[SandboxCase] = []
        for gene_dir in sorted(genes_dir.iterdir()):
            if not gene_dir.is_dir():
                continue
            for case_dir in sorted(gene_dir.iterdir()):
                manifest = case_dir / "manifest.yaml"
                if not manifest.exists():
                    continue
                data = yaml.safe_load(manifest.read_text(encoding="utf-8"))
                cases.append(self._case_from_manifest(case_dir, data))
        return cases

    def _case_from_manifest(self, case_dir: Path, data: dict) -> SandboxCase:
        gene = data.get("gene", case_dir.parent.name)
        spdi = data.get("spdi", "")
        if not spdi:
            raise SandboxError(f"manifest without spdi: {case_dir}")
        case_id = f"{gene}:{spdi}"

        bams: dict[str, Path] = {}
        bams_dir = case_dir / "bams"
        if bams_dir.exists():
            for sample in ("father", "mother", "proband"):
                p = bams_dir / f"{sample}.bam"
                if p.exists():
                    bams[sample] = p

        gt = data.get("ground_truth", {}) or {}
        chrom = gt.get("chrom", "")
        zygosity = gt.get("zygosity", {}) or {}

        return SandboxCase(
            case_id=case_id,
            gene=gene,
            spdi=spdi,
            title=data.get("title", ""),
            clinical_text=data.get("clinical_text", ""),
            chrom=chrom,
            zygosity=zygosity,
            hpo_expected=data.get("hpo_expected", []) or [],
            path=case_dir.resolve(),
            bams=bams,
        )

    # --- conversion to PGA internal types --------------------------------

    def to_case_manifest(self, case: SandboxCase, *, roi_padding: int = 5000) -> CaseManifest:
        """Build a CaseManifest from a SandboxCase.

        The pedigree is inferred from the trio structure
        (father / mother / proband) plus the zygosity map.
        Affected status: proband if its zygosity is HOM_DEL or HOM_ALT,
        parents unaffected otherwise.
        """
        # Pedigree
        parents_present = "father" in case.bams and "mother" in case.bams
        proband_affected = case.zygosity.get("proband") in ("HOM_DEL", "HOM_ALT")

        pedigree: list[IndividualSpec] = []
        if parents_present:
            pedigree.append(IndividualSpec(
                sample="father", sex="M", affected=False,
            ))
            pedigree.append(IndividualSpec(
                sample="mother", sex="F", affected=False,
            ))
            pedigree.append(IndividualSpec(
                sample="proband", sex="U", affected=proband_affected,
                father="father", mother="mother",
            ))
            proband = "proband"
        else:
            # Single sample; treat as a lone proband
            sample = next(iter(case.bams), "proband")
            pedigree.append(IndividualSpec(
                sample=sample, sex="U", affected=True,
            ))
            proband = sample

        # ROI from SPDI + chrom
        rois: list[GenomicRange] = []
        if case.chrom and case.spdi:
            try:
                _refseq, pos, del_len, ins_len_seq = _parse_spdi(case.spdi)
                start = max(1, pos - roi_padding)
                end = pos + max(del_len, len(ins_len_seq)) + roi_padding
                rois.append(GenomicRange(
                    chrom=case.chrom,
                    start=start,
                    end=end,
                    label=case.gene,
                ))
            except SandboxError:
                pass  # no ROI if SPDI is malformed

        # HPO terms
        hpo_terms: dict[str, list[str]] = {}
        if case.hpo_expected:
            hpo_terms[proband] = [h.get("code", "") for h in case.hpo_expected if h.get("code")]

        return CaseManifest(
            case_id=case.case_id.replace(":", "_"),
            source="clinical",
            pedigree=pedigree,
            proband=proband,
            affected_samples=[proband] if proband_affected else [],
            candidate_rois=rois,
            hpo_terms=hpo_terms,
            phenotype_text={proband: case.clinical_text} if case.clinical_text else {},
            candidate_genes=[case.gene] if case.gene else [],
            consanguinity=False,
            confidence=1.0,
            extractor="sandbox",
            confirmed_by=None,
        )

