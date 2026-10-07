"""Fixtures compartidas para los tests de integración.

Un test de integración asume que los entornos micromamba
(pga-core, pga-hts, etc.) están creados y disponibles en
~/.pga/envs/. No intenta crearlos: si faltan, el test falla
con un mensaje claro. Los tests unitarios no dependen de esto.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
MAKE_DEMO_TRIO = REPO_ROOT / "scripts" / "make_demo_trio.py"


@pytest.fixture(scope="function")
def demo_trio(tmp_path: Path) -> dict:
    """Genera un trío sintético (father/mother/proband) en tmp_path.

    Devuelve un dict con:
      - out_dir: directorio con todos los ficheros generados
      - father_bam, mother_bam, proband_bam: rutas a los BAMs
      - ref_fa: ruta a la referencia FASTA
      - case_json: ruta al case.json (CaseManifest)
      - case: dict cargado del case.json
      - region: "chr1:1000-1200" (ROI del trío)

    El trío tiene un variant recesivo en chr1:1100 (A>G):
      - father: 0/1
      - mother: 0/1
      - proband: 1/1
    """
    if not MAKE_DEMO_TRIO.exists():
        pytest.skip(f"make_demo_trio.py no encontrado en {MAKE_DEMO_TRIO}")

    out_dir = tmp_path / "demo_trio"

    result = subprocess.run(
        [sys.executable, str(MAKE_DEMO_TRIO), str(out_dir)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"make_demo_trio.py falló (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )

    case_json = out_dir / "case.json"
    if not case_json.exists():
        pytest.fail(f"case.json no generado en {out_dir}")

    return {
        "out_dir": out_dir,
        "father_bam": out_dir / "father.bam",
        "mother_bam": out_dir / "mother.bam",
        "proband_bam": out_dir / "proband.bam",
        "ref_fa": out_dir / "ref.fa",
        "case_json": case_json,
        "case": json.loads(case_json.read_text()),
        "region": "chr1:1000-1200",
        "variant": {"chrom": "chr1", "pos": 1100, "ref": "A", "alt": "G"},
    }

@pytest.fixture(scope="function")
def synthetic_sandbox(tmp_path: Path) -> dict:
    """Build a minimal Sandbox case in tmp_path with 3 real BAMs.

    Creates a SNV trio: father/mother HET, proband HOM_ALT, at
    chr1:500 (A>G) on a synthetic reference of 1000 bp.

    Returns paths to the case directory, the reference, and the
    expected genotypes.
    """
    import subprocess
    import pysam

    # --- 1. Reference: 1000 bp of A, chr1 ---
    case_dir = tmp_path / "genes" / "SYNTH" / "NC_000001.11_500_A_G"
    bams_dir = case_dir / "bams"
    bams_dir.mkdir(parents=True, exist_ok=True)

    ref_fa = case_dir / "ref.fa"
    seq = "A" * 1000
    with ref_fa.open("w") as f:
        f.write(">chr1\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
    pysam.faidx(str(ref_fa))

    # --- 2. Manifest ---
    import yaml
    manifest = {
        "gene": "SYNTH",
        "spdi": "NC_000001.11:500:A:G",
        "title": "Synthetic SNV trio",
        "clinical_text": "Synthetic test case.",
        "ground_truth": {
            "chrom": "chr1",
            "zygosity": {
                "father": "HET",
                "mother": "HET",
                "proband": "HOM_ALT",
            },
        },
        "hpo_expected": [],
    }
    (case_dir / "manifest.yaml").write_text(yaml.safe_dump(manifest))

    # --- 3. Three BAMs ---
    def _make_bam(sample: str, gt: str) -> Path:
        """Generate a BAM with 50 reads covering chr1:450-550."""
        bam_path = bams_dir / f"{sample}.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
            "RG": [{"ID": "rg1", "SM": sample, "PL": "ILLUMINA"}],
        }
        import random
        rng = random.Random(42)
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
            for i in range(100):
                start = rng.randint(400, 500)
                read_len = 100
                # Build the read from the reference
                bases = ["A"] * read_len
                # Apply variant at chr1:500 if this read carries it
                pos_in_read = 500 - start  # 0-based offset of pos 500
                if 0 <= pos_in_read < read_len:
                    if gt == "0/1" and rng.random() < 0.5:
                        bases[pos_in_read] = "G"
                    elif gt == "1/1":
                        bases[pos_in_read] = "G"
                a = pysam.AlignedSegment()
                a.query_name = f"r{i:05d}"
                a.query_sequence = "".join(bases)
                a.flag = 0
                a.reference_id = 0
                a.reference_start = start
                a.mapping_quality = 60
                a.cigar = ((0, read_len),)
                a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
                a.set_tag("RG", "rg1")
                out.write(a)

        # Sort + index
        sorted_bam = bams_dir / f"{sample}.sorted.bam"
        pysam.sort("-o", str(sorted_bam), str(bam_path))
        sorted_bam.rename(bam_path)
        pysam.index(str(bam_path))
        return bam_path

    _make_bam("father", "0/1")
    _make_bam("mother", "0/1")
    _make_bam("proband", "1/1")

    return {
        "case_dir": case_dir,
        "ref_fa": ref_fa,
        "expected_genotypes": {
            "father": "0/1",
            "mother": "0/1",
            "proband": "1/1",
        },
    }

