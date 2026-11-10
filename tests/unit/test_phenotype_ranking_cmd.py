"""LIRICAL: formato numérico fijo y modo fenotipo por defecto."""
from __future__ import annotations

from pathlib import Path

from publicgenomicagent.tools.base import PhenotypeRankingInput
from publicgenomicagent.tools.phenotype_ranking import _build_command


def test_lirical_command_forces_english_numbers_and_phenotype_only(tmp_path):
    inp = PhenotypeRankingInput(hpo_ids=["HP:0000090", "HP:0000083"], output_dir=tmp_path)
    cmd = [str(c) for c in _build_command(inp, Path("/x/lirical-cli-2.4.1.jar"), Path("/x/data"), "p")]
    assert cmd.index("-Duser.language=en") < cmd.index("-jar")
    assert "-Duser.country=US" in cmd and "--vcf" not in cmd
    assert cmd[cmd.index("-p") + 1] == "HP:0000090,HP:0000083"
