from __future__ import annotations

import gzip
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from ..env.micromamba import bin_path
from ..env.runtime import ToolRuntime
from .base import (
    IndividualSpec,
    MendelianFilterInput,
    MendelianFilterOutput,
    PlinkValidateInput,
    PlinkValidateOutput,
)


# ---------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------

def _bcftools() -> str:
    p = bin_path("pga-hts", "bcftools")
    if not p.exists():
        raise FileNotFoundError(f"bcftools no encontrado: {p}")
    return str(p)


def _plink() -> str:
    p = bin_path("pga-mendelian", "plink")
    if not p.exists():
        raise FileNotFoundError(f"plink no encontrado: {p}")
    return str(p)


@dataclass
class SampleData:
    name: str
    gt: str          # "0/1", "1/1", "0/0", "./.", etc.
    dp: int | None


def _parse_gt(gt_field: str) -> tuple[str, int | None]:
    """Extrae GT y DP de un campo FORMAT/muestra."""
    if ":" not in gt_field:
        return gt_field, None
    parts = gt_field.split(":")
    gt = parts[0]
    # No tenemos acceso a los índices del FORMAT aquí; se resuelve en el caller
    return gt, None


def _sample_gt(fields: list[str], samples: list[str], sample: str,
               fmt_keys: list[str]) -> tuple[str, int | None]:
    """Devuelve (GT, DP) para un sample concreto."""
    idx = samples.index(sample)
    sample_field = fields[9 + idx].split(":")
    gt = sample_field[0] if sample_field else "./."
    dp = None
    if "DP" in fmt_keys:
        dp_idx = fmt_keys.index("DP")
        if dp_idx < len(sample_field) and sample_field[dp_idx].isdigit():
            dp = int(sample_field[dp_idx])
    return gt, dp


def _iter_vcf(vcf: Path):
    """Iterador sobre líneas de un VCF, gz o plano."""
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt") as f:
        for line in f:
            yield line


def _write_vcf_header(out_f, header_lines: list[str]) -> None:
    for line in header_lines:
        if line.startswith("##"):
            out_f.write(line)
    for line in header_lines:
        if line.startswith("#CHROM"):
            out_f.write(line)
            break


# ---------------------------------------------------------------------
# Escritura de sub-VCFs
# ---------------------------------------------------------------------

def _write_subvcf(
    out_path: Path,
    header_lines: list[str],
    matching_lines: list[str],
) -> None:
    """Escribe un VCF con el header del original y las líneas filtradas.

    Escribe primero un VCF plano, luego comprime con `bcftools view -Oz`
    (que usa BGZF, necesario para indexar con tabix). Si no hay
    variantes, deja el VCF vacío sin indexar.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plain = out_path.with_suffix(".plain.vcf")
    with plain.open("w") as f:
        _write_vcf_header(f, header_lines)
        for line in matching_lines:
            f.write(line)

    # Comprimir a BGZF con bcftools
    subprocess.run(
        [_bcftools(), "view", str(plain), "-Oz", "-o", str(out_path)],
        check=True, capture_output=True,
    )
    plain.unlink()

    # Indexar solo si hay variantes
    if matching_lines:
        subprocess.run(
            [_bcftools(), "index", "-t", str(out_path)],
            check=True, capture_output=True,
        )


# ---------------------------------------------------------------------
# Reglas de segregación
# ---------------------------------------------------------------------

def _parse_fmt(fmt_field: str) -> list[str]:
    """Divide el campo FORMAT en claves."""
    return fmt_field.split(":")


def _is_het(gt: str) -> bool:
    return gt in ("0/1", "1/0", "0|1", "1|0")


def _is_hom_alt(gt: str) -> bool:
    return gt in ("1/1", "1|1")


def _is_hom_ref(gt: str) -> bool:
    return gt in ("0/0", "0|0")


def _is_alt_carrier(gt: str) -> bool:
    """Portador: het o hom alt."""
    return _is_het(gt) or _is_hom_alt(gt)


def _passes_quality(
    fields: list[str], fmt_keys: list[str],
    samples: list[str], proband: str,
    min_dp: int, min_qual: int,
) -> bool:
    """Verifica que QUAL y el DP del probando cumplen los mínimos.

    El DP se exige solo al probando. Los padres pueden tener DP=0
    o ./., por ejemplo si no portan la variante.
    """
    # QUAL
    try:
        qual = float(fields[5]) if fields[5] != "." else 0.0
    except ValueError:
        qual = 0.0
    if qual < min_qual:
        return False

    # DP del probando
    if "DP" not in fmt_keys:
        return True
    if proband not in samples:
        return True
    dp_idx = fmt_keys.index("DP")
    idx = samples.index(proband)
    sample_field = fields[9 + idx].split(":")
    if dp_idx < len(sample_field) and sample_field[dp_idx].isdigit():
        if int(sample_field[dp_idx]) < min_dp:
            return False
    return True


def _filter_de_novo(
    vcf_lines: list[str], header_lines: list[str],
    samples: list[str], child: str, father: str, mother: str,
    min_dp: int, min_qual: int,
) -> tuple[list[str], int]:
    """Padre=0/0 AND madre=0/0 AND hijo portador."""
    matches: list[str] = []
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        fmt_keys = _parse_fmt(fields[8])
        if not _passes_quality(fields, fmt_keys, samples, child, min_dp, min_qual):
            continue
        gt_child, _ = _sample_gt(fields, samples, child, fmt_keys)
        gt_father, _ = _sample_gt(fields, samples, father, fmt_keys)
        gt_mother, _ = _sample_gt(fields, samples, mother, fmt_keys)
        if (_is_hom_ref(gt_father) and _is_hom_ref(gt_mother)
                and _is_alt_carrier(gt_child)):
            matches.append(line)
    return matches, len(matches)


def _filter_auto_rec_hom(
    vcf_lines: list[str], header_lines: list[str],
    samples: list[str], child: str, father: str, mother: str,
    min_dp: int, min_qual: int,
) -> tuple[list[str], int]:
    """Padre=0/1 AND madre=0/1 AND hijo=1/1."""
    matches: list[str] = []
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        fmt_keys = _parse_fmt(fields[8])
        if not _passes_quality(fields, fmt_keys, samples, child, min_dp, min_qual):
            continue
        gt_child, _ = _sample_gt(fields, samples, child, fmt_keys)
        gt_father, _ = _sample_gt(fields, samples, father, fmt_keys)
        gt_mother, _ = _sample_gt(fields, samples, mother, fmt_keys)
        if (_is_het(gt_father) and _is_het(gt_mother)
                and _is_hom_alt(gt_child)):
            matches.append(line)
    return matches, len(matches)


def _filter_auto_dom(
    vcf_lines: list[str], header_lines: list[str],
    samples: list[str], child: str, father: str, mother: str,
    affected_samples: list[str],
    min_dp: int, min_qual: int,
) -> tuple[list[str], int]:
    """Hijo portador y al menos un progenitor portador, con al menos un
    progenitor afectado portador."""
    matches: list[str] = []
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        fmt_keys = _parse_fmt(fields[8])
        if not _passes_quality(fields, fmt_keys, samples, child, min_dp, min_qual):
            continue
        gt_child, _ = _sample_gt(fields, samples, child, fmt_keys)
        gt_father, _ = _sample_gt(fields, samples, father, fmt_keys)
        gt_mother, _ = _sample_gt(fields, samples, mother, fmt_keys)
        if not _is_alt_carrier(gt_child):
            continue
        # Al menos un progenitor portador
        parent_carriers: list[str] = []
        if _is_alt_carrier(gt_father):
            parent_carriers.append(father)
        if _is_alt_carrier(gt_mother):
            parent_carriers.append(mother)
        if not parent_carriers:
            continue
        # Si hay lista de afectados, exigir que al menos un portador esté
        # en affected_samples. Si está vacía, aceptar cualquier portador.
        if affected_samples:
            if not any(p in affected_samples for p in parent_carriers):
                continue
        matches.append(line)
    return matches, len(matches)


def _filter_x_linked_rec(
    vcf_lines: list[str], header_lines: list[str],
    samples: list[str], child: str, father: str, mother: str,
    child_sex: str, min_dp: int, min_qual: int,
) -> tuple[list[str], int]:
    """Hijo varón hemicigoto (1/1), madre portadora (0/1), padre 0/0.

    Solo se aplica a variantes en el cromosoma X. Si no hay cromosoma X
    en el VCF, devuelve lista vacía.
    """
    matches: list[str] = []
    if child_sex != "M":
        return matches, 0
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        chrom = fields[0].lower().replace("chr", "")
        if chrom not in ("x", "23"):
            continue
        fmt_keys = _parse_fmt(fields[8])
        if not _passes_quality(fields, fmt_keys, samples, child, min_dp, min_qual):
            continue
        gt_child, _ = _sample_gt(fields, samples, child, fmt_keys)
        gt_father, _ = _sample_gt(fields, samples, father, fmt_keys)
        gt_mother, _ = _sample_gt(fields, samples, mother, fmt_keys)
        if (_is_hom_alt(gt_child) and _is_het(gt_mother)
                and _is_hom_ref(gt_father)):
            matches.append(line)
    return matches, len(matches)


# ---------------------------------------------------------------------
# Entry point: mendelian_filter
# ---------------------------------------------------------------------

def _load_vcf(vcf: Path) -> tuple[list[str], list[str], list[str]]:
    """Devuelve (header_lines, data_lines, samples)."""
    header: list[str] = []
    data: list[str] = []
    samples: list[str] = []
    for line in _iter_vcf(vcf):
        if line.startswith("##"):
            header.append(line)
        elif line.startswith("#CHROM"):
            header.append(line)
            fields = line.rstrip("\n").split("\t")
            samples = fields[9:]
        elif line.strip() and not line.startswith("#"):
            data.append(line)
    if not samples:
        raise ValueError(f"VCF sin samples: {vcf}")
    return header, data, samples


def mendelian_filter(inp: MendelianFilterInput) -> MendelianFilterOutput:
    """Aplica reglas de segregación mendeliana sobre un VCF de trío.

    Modelos implementados:
    - de_novo
    - auto_rec_hom (recesivo autosómico homocigoto)
    - auto_dom (dominante autosómico)
    - x_linked_rec (ligado a X recesivo, solo si hay cromosoma X)
    """
    vcf = Path(inp.trio_vcf).expanduser().resolve()
    outdir = Path(inp.output_dir).expanduser().resolve()

    if not vcf.exists():
        raise FileNotFoundError(f"VCF no encontrado: {vcf}")
    outdir.mkdir(parents=True, exist_ok=True)

    # Extraer padre y madre del probando
    father, mother = None, None
    for ind in inp.pedigree:
        if ind.sample == inp.proband:
            father = ind.father
            mother = ind.mother
            break

    if father is None or mother is None:
        raise ValueError(
            f"El probando '{inp.proband}' no tiene padre y madre definidos "
            f"en el pedigrí. Un análisis de trío requiere ambos."
        )

    # Cargar el VCF
    header, data, samples = _load_vcf(vcf)

    # Verificar que los tres samples están
    for s in (inp.proband, father, mother):
        if s not in samples:
            raise ValueError(
                f"Sample '{s}' no está en el VCF. Disponibles: {samples}"
            )

    # Derivar affected_samples del pedigrí.
    # Si filter_by_affected es False, pasar lista vacía a auto_dom
    # para relajar el modelo (penetrancia incompleta, variantes de
    # riesgo).
    if inp.filter_by_affected:
        affected_samples = [ind.sample for ind in inp.pedigree if ind.affected]
    else:
        affected_samples = []

    # Aplicar los filtros
    de_novo_lines, n_dn = _filter_de_novo(
        data, header, samples, inp.proband, father, mother,
        inp.min_dp, inp.min_qual,
    )
    rec_lines, n_rec = _filter_auto_rec_hom(
        data, header, samples, inp.proband, father, mother,
        inp.min_dp, inp.min_qual,
    )
    dom_lines, n_dom = _filter_auto_dom(
        data, header, samples, inp.proband, father, mother,
        affected_samples, inp.min_dp, inp.min_qual,
    )

    # Sexo del probando para X-linked
    child_sex = "U"
    for ind in inp.pedigree:
        if ind.sample == inp.proband:
            child_sex = ind.sex
            break
    xrec_lines, n_xrec = _filter_x_linked_rec(
        data, header, samples, inp.proband, father, mother,
        child_sex, inp.min_dp, inp.min_qual,
    )

    # Escribir los sub-VCFs
    de_novo_vcf = outdir / "de_novo.vcf.gz"
    rec_vcf = outdir / "auto_rec_hom.vcf.gz"
    dom_vcf = outdir / "auto_dom.vcf.gz"
    xrec_vcf = outdir / "x_linked_rec.vcf.gz"

    _write_subvcf(de_novo_vcf, header, de_novo_lines)
    _write_subvcf(rec_vcf, header, rec_lines)
    _write_subvcf(dom_vcf, header, dom_lines)
    _write_subvcf(xrec_vcf, header, xrec_lines)

    # Report
    report = {
        "proband": inp.proband,
        "father": father,
        "mother": mother,
        "child_sex": child_sex,
        "affected_samples": affected_samples,
        "counts": {
            "de_novo": n_dn,
            "auto_rec_hom": n_rec,
            "auto_dom": n_dom,
            "x_linked_rec": n_xrec,
        },
        "total_input_variants": len(data),
        "min_dp": inp.min_dp,
        "min_qual": inp.min_qual,
    }
    report_path = outdir / "mendelian_report.json"
    report_path.write_text(json.dumps(report, indent=2))

    return MendelianFilterOutput(
        ok=True,
        tool="mendelian_filter",
        message=(
            f"de_novo={n_dn}  auto_rec_hom={n_rec}  "
            f"auto_dom={n_dom}  x_linked_rec={n_xrec}"
        ),
        outputs={
            "de_novo": str(de_novo_vcf),
            "auto_rec_hom": str(rec_vcf),
            "auto_dom": str(dom_vcf),
            "x_linked_rec": str(xrec_vcf),
            "report": str(report_path),
        },
        de_novo_vcf=de_novo_vcf,
        auto_rec_hom_vcf=rec_vcf,
        auto_dom_vcf=dom_vcf,
        x_linked_rec_vcf=xrec_vcf,
        report_json=report_path,
        counts=report["counts"],
    )


# ---------------------------------------------------------------------
# Entry point: plink_validate
# ---------------------------------------------------------------------

def _write_fam(pedigree: list[IndividualSpec], out_fam: Path) -> None:
    """Escribe un .fam desde el pedigrí.

    Formato .fam: FID IID PID MID sex phenotype
    - sex: 1=M, 2=F, 0=U
    - phenotype: 1=unaffected, 2=affected, 0=missing
    """
    sex_map = {"M": "1", "F": "2", "U": "0"}
    with out_fam.open("w") as f:
        for ind in pedigree:
            fid = ind.sample
            iid = ind.sample
            pid = ind.father or "0"
            mid = ind.mother or "0"
            sex = sex_map.get(ind.sex, "0")
            pheno = "2" if ind.affected else "1"
            f.write(f"{fid}\t{iid}\t{pid}\t{mid}\t{sex}\t{pheno}\n")


def _vcf_to_plink(vcf: Path, out_prefix: Path) -> None:
    """Convierte VCF a formato PLINK (bed/bim/fam).

    Usa --allow-extra-chr para aceptar nombres de cromosoma no
    estándar (por ejemplo, chr2_roi, o subregiones locales). Usa
    --double-id para aceptar sample IDs sin formato FID/IID.
    """
    subprocess.run(
        [_plink(), "--vcf", str(vcf),
         "--make-bed",
         "--out", str(out_prefix),
         "--allow-no-sex",
         "--allow-extra-chr",
         "--double-id",
         "--silent"],
        check=True, capture_output=True,
    )


def plink_validate(inp: PlinkValidateInput) -> PlinkValidateOutput:
    """Valida consistencia mendeliana con PLINK.

    - Convierte VCF a PLINK.
    - Sobreescribe el .fam con el pedigrí del manifiesto.
    - Corre --mendel para detectar errores mendelianos.
    - Corre --genome para calcular IBD entre individuos.
    """
    vcf = Path(inp.trio_vcf).expanduser().resolve()
    outdir = Path(inp.output_dir).expanduser().resolve()

    if not vcf.exists():
        raise FileNotFoundError(f"VCF no encontrado: {vcf}")
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = outdir / "plink_trio"

    # Convertir VCF a PLINK
    _vcf_to_plink(vcf, prefix)

    # Sobreescribir .fam con el pedigrí del manifiesto
    _write_fam(inp.pedigree, Path(str(prefix) + ".fam"))

    # --mendel: errores mendelianos
    mendel_path = outdir / "plink_mendel.txt"
    with mendel_path.open("w") as f:
        subprocess.run(
            [_plink(), "--bfile", str(prefix),
             "--mendel",
             "--out", str(outdir / "plink_mendel_run"),
             "--allow-no-sex",
             "--allow-extra-chr",
             "--silent"],
            check=True, stdout=f, stderr=subprocess.STDOUT,
        )

    # Contar errores mendelianos
    n_errors = 0
    mendel_summary = outdir / "plink_mendel_run.mendel"
    if mendel_summary.exists():
        for line in mendel_summary.read_text().splitlines():
            if "NONE" in line.upper() or not line.strip():
                continue
            n_errors += 1
    else:
        # PLINK no crea el .mendel si no hay errores.
        # Creamos un placeholder para dejar constancia de que se ejecutó.
        mendel_summary.write_text(
            "# PLINK no reportó errores mendelianos (fichero vacío).\n"
        )

    # --genome: IBD entre individuos
    ibd_path = outdir / "plink_ibd.genome"
    subprocess.run(
        [_plink(), "--bfile", str(prefix),
         "--genome",
         "--out", str(outdir / "plink_ibd_run"),
         "--allow-no-sex",
         "--allow-extra-chr",
         "--silent"],
        check=True, capture_output=True,
    )
    ibd_run = outdir / "plink_ibd_run.genome"
    if ibd_run.exists():
        ibd_run.replace(ibd_path)

    # Report
    report = {
        "proband": inp.proband,
        "n_samples": len(inp.pedigree),
        "n_mendel_errors": n_errors,
        "bed_prefix": str(prefix),
        "mendel_file": str(mendel_summary),
        "ibd_file": str(ibd_path),
    }
    report_path = outdir / "plink_validate_report.json"
    report_path.write_text(json.dumps(report, indent=2))

    return PlinkValidateOutput(
        ok=True,
        tool="plink_validate",
        message=(
            f"PLINK: {len(inp.pedigree)} samples, "
            f"{n_errors} errores mendelianos"
        ),
        outputs={
            "bed": str(prefix),
            "mendel": str(mendel_summary),
            "ibd": str(ibd_path),
            "report": str(report_path),
        },
        bed_prefix=prefix,
        mendel_errors=mendel_summary,
        ibd_report=ibd_path,
        report_json=report_path,
        n_mendel_errors=n_errors,
    )
