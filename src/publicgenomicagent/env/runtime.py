from __future__ import annotations

import re
import subprocess

from .micromamba import bin_path, env_exists
from .registry import Registry, ToolSpec

_VERSION_RE = re.compile(r"\d+\.\d+(?:\.\d+)?")


class ToolRuntime:
    """Resuelve y ejecuta herramientas a través de sus entornos aislados.

    Una vez creado un entorno, sus binarios se invocan directamente
    (no vía `micromamba run`), lo cual evita el bug de micromamba 2.x
    que intercepta `--version` y permite lanzar comandos más rápido.
    """

    def __init__(self, registry: Registry):
        self.registry = registry

    def spec(self, tool: str) -> ToolSpec:
        if tool not in self.registry.tools:
            raise KeyError(f"Herramienta '{tool}' no registrada")
        return self.registry.tools[tool]

    def ensure_env_for(self, tool: str) -> None:
        spec = self.spec(tool)
        if not env_exists(spec.env):
            raise RuntimeError(
                f"El entorno '{spec.env}' no está creado. "
                f"Ejecuta: pga env bootstrap {spec.env}"
            )

    def run(self, tool: str, args: list[str], **kwargs) -> subprocess.CompletedProcess:
        spec = self.spec(tool)
        self.ensure_env_for(tool)
        binary = bin_path(spec.env, spec.binary)
        if not binary.exists():
            raise FileNotFoundError(f"Binario no encontrado: {binary}")
        return subprocess.run([str(binary), *args], check=True, **kwargs)

    def capture(self, tool: str, args: list[str]) -> str:
        result = self.run(tool, args, capture_output=True, text=True)
        # Unificamos stdout y stderr: samtools/bcftools escriben --version en stderr.
        return ((result.stdout or "") + "\n" + (result.stderr or "")).strip()

    def version(self, tool: str) -> str:
        spec = self.spec(tool)
        raw = self.capture(tool, [spec.version_flag])
        for line in raw.splitlines():
            m = _VERSION_RE.search(line)
            if m:
                return m.group(0)
        return raw.splitlines()[0] if raw else "desconocida"
