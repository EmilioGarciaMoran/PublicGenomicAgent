from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class EnvSpec:
    name: str
    manifest: Path
    description: str
    installs_self: bool = False


@dataclass(frozen=True)
class ToolSpec:
    name: str
    env: str
    binary: str
    version_flag: str
    min_version: str | None


@dataclass(frozen=True)
class Registry:
    envs: dict[str, EnvSpec]
    tools: dict[str, ToolSpec]


def load_registry(registry_path: Path) -> Registry:
    if not registry_path.exists():
        raise FileNotFoundError(f"registry.yaml no encontrado en {registry_path}")

    raw = yaml.safe_load(registry_path.read_text())

    envs = {
        name: EnvSpec(
            name=name,
            manifest=Path(spec["manifest"]),
            description=spec.get("description", ""),
            installs_self=spec.get("installs_self", False),
        )
        for name, spec in raw.get("envs", {}).items()
    }

    tools = {
        name: ToolSpec(
            name=name,
            env=spec["env"],
            binary=spec["binary"],
            version_flag=spec.get("version_flag", "--version"),
            min_version=spec.get("min_version"),
        )
        for name, spec in raw.get("tools", {}).items()
    }

    return Registry(envs=envs, tools=tools)
