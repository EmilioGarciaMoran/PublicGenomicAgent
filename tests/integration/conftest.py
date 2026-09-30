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
