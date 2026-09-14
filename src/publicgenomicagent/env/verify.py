from __future__ import annotations

from dataclasses import dataclass

from .manifests import read_lock, sha256_of_file
from .micromamba import env_exists
from .registry import Registry
from .runtime import ToolRuntime


@dataclass
class CheckResult:
    tool: str
    env: str
    status: str  # ok | missing_env | missing_binary | version_mismatch | stale_manifest
    detail: str = ""


def verify_all(registry: Registry) -> list[CheckResult]:
    runtime = ToolRuntime(registry)
    results: list[CheckResult] = []

    for env_name, spec in registry.envs.items():
        if not env_exists(env_name):
            results.append(CheckResult(
                tool="-", env=env_name, status="missing_env",
                detail="entorno no creado",
            ))
            continue

        lock = read_lock(env_name)
        current_hash = sha256_of_file(spec.manifest)
        if lock is None or lock.get("manifest_hash") != current_hash:
            results.append(CheckResult(
                tool="-", env=env_name, status="stale_manifest",
                detail="manifiesto cambió desde el último bootstrap",
            ))

    for tool_name, tool_spec in registry.tools.items():
        try:
            version = runtime.version(tool_name)
            results.append(CheckResult(
                tool=tool_name, env=tool_spec.env, status="ok",
                detail=version,
            ))
        except RuntimeError as e:
            results.append(CheckResult(
                tool=tool_name, env=tool_spec.env,
                status="missing_env", detail=str(e),
            ))
        except Exception as e:  # noqa: BLE001
            results.append(CheckResult(
                tool=tool_name, env=tool_spec.env,
                status="missing_binary", detail=str(e),
            ))

    return results
