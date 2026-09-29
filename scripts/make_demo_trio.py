#!/usr/bin/env python3
"""Genera un trío sintético con un variant recesivo plantado.

Caso: enfermedad autosómica recesiva.
  - Father: het (0/1)
  - Mother: het (0/1)
  - Proband: hom alt (1/1), afectado

Locus: chr1:1100, A>G, en un ROI chr1:1000-1200 (label "DEMO1").
Referencia: tiny (chr1:1-10000 de puras A) generada al vuelo.

Uso:
    python scripts/make_demo_trio.py [<out_dir>]
    # por defecto: results/demo_trio
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures"

# Coordenadas del variant
CHROM = "chr1"
POS = 1100
REF = "A"
ALT = "G"

# Genotipos mendelianos del trío
TRIO = [
    {"sample": "father",  "sex": "M", "affected": False, "gt": "0/1"},
    {"sample": "mother",  "sex": "F", "affected": False, "gt": "0/1"},
    {"sample": "proband", "sex": "U", "affected": True,  "gt": "1/1",
     "father": "father", "mother": "mother"},
]

# ROI
ROI = {"chrom": CHROM, "start": 1000, "end": 1200, "label": "DEMO1"}

# Cobertura de los BAMs sintéticos
COVERAGE = 60


def run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise SystemExit(f"ERROR: comando falló con código {result.returncode}")


def main() -> int:
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / "results" / "demo_trio"
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_fa = out_dir / "ref.fa"

    # 1. Generar referencia tiny
    print(f"[1/4] Generando referencia tiny en {ref_fa}")
    run([
        sys.executable,
        str(FIXTURES / "make_tiny_ref.py"),
        str(ref_fa),
    ])

    # 2. Generar BAMs del trío
    print(f"[2/4] Generando 3 BAMs (cobertura {COVERAGE}x)")
    for ind in TRIO:
        bam = out_dir / f"{ind['sample']}.bam"
        print(f"      {ind['sample']} ({ind['gt']}) → {bam.name}")
        run([
            sys.executable,
            str(FIXTURES / "make_bam_with_variant.py"),
            str(ref_fa),
            CHROM,
            str(POS),
            REF,
            ALT,
            ind["sample"],
            str(bam),
            str(COVERAGE),
            ind["gt"],
        ])

    # 3. Escribir case.json
    print(f"[3/4] Escribiendo case.json")
    case = {
        "case_id": "DEMO_TRIO",
        "source": "clinical",
        "pedigree": [
            {
                "sample": ind["sample"],
                "sex": ind["sex"],
                "affected": ind["affected"],
                "father": ind.get("father"),
                "mother": ind.get("mother"),
            }
            for ind in TRIO
        ],
        "proband": "proband",
        "affected_samples": ["proband"],
        "candidate_rois": [ROI],
        "hpo_terms": {},
        "phenotype_text": {},
        "candidate_genes": [],
        "consanguinity": True,
        "confidence": 1.0,
        "extractor": "manual",
        "confirmed_by": None,
    }
    case_path = out_dir / "case.json"
    case_path.write_text(json.dumps(case, indent=2), encoding="utf-8")

    # 4. Resumen
    print(f"[4/4] OK. Ficheros en {out_dir}:")
    for f in sorted(out_dir.iterdir()):
        print(f"      {f.name}  ({f.stat().st_size} bytes)")

    print()
    print("Siguiente paso:")
    print(f"  pga run-trio \\")
    print(f"    --father {out_dir}/father.bam \\")
    print(f"    --mother {out_dir}/mother.bam \\")
    print(f"    --proband {out_dir}/proband.bam \\")
    print(f"    --reference {ref_fa} \\")
    print(f"    --region {CHROM}:{ROI['start']}-{ROI['end']} \\")
    print(f"    --case-id DEMO_TRIO --label DEMO1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
