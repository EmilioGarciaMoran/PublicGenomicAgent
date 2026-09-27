from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

from ..env.micromamba import bin_path
from .base import ExtractHPOInput, ExtractHPOOutput


# ---------------------------------------------------------------------
# Configuración
# ---------------------------------------------------------------------

DEFAULT_RDMA_CACHE = Path.home() / ".pga" / "cache" / "rdma"
DEFAULT_MODEL_CACHE = Path.home() / ".pga" / "cache" / "models"

HPO_EMBEDDINGS = "vector_stores/G2GHPO_metadata_medembed.npy"
LAB_EMBEDDINGS = "tools/lab_tables_medembed_sm.npy"


def _python() -> str:
    """Python del entorno pga-phenotype."""
    p = bin_path("pga-phenotype", "python")
    if not p.exists():
        raise FileNotFoundError(
            f"Python de pga-phenotype no encontrado: {p}. "
            "Ejecuta: pga env bootstrap pga-phenotype"
        )
    return str(p)


def _resolve_cache_dir(inp: ExtractHPOInput) -> Path:
    c = inp.rdma_cache_dir or DEFAULT_RDMA_CACHE
    c = Path(c).expanduser().resolve()
    if not (c / HPO_EMBEDDINGS).exists():
        raise FileNotFoundError(
            f"Vector stores de RDMA no encontrados en {c}. "
            f"Esperado: {c / HPO_EMBEDDINGS}"
        )
    return c


def _resolve_model_cache(inp: ExtractHPOInput) -> Path:
    c = inp.model_cache_dir or DEFAULT_MODEL_CACHE
    c = Path(c).expanduser().resolve()
    c.mkdir(parents=True, exist_ok=True)
    return c


def _text_hash(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


# ---------------------------------------------------------------------
# Runner: ejecuta RDMA en subproceso para aislar el LLM
# ---------------------------------------------------------------------

_RUNNER = r'''
import json
import sys
from pathlib import Path

# Parámetros vía argv
cfg = json.loads(sys.argv[1])

from rdma.hpo.extractor import PhenotypeExtractor
from rdma.hpo.verifier import HPOVerifier
from rdma.hpo.matcher import HPOMatcher

# Verificación previa: backend local requiere GPU
backend = cfg["backend"]
if backend == "local":
    import torch
    device = cfg["device"]
    if device == "auto":
        if not torch.cuda.is_available():
            raise RuntimeError(
                "Backend 'local' requiere GPU NVIDIA con CUDA. "
                "No se detectó CUDA disponible. "
                "Alternativas: --backend openrouter o --backend api, "
                "o usar una máquina con GPU (>=20 GB VRAM para Mistral 24B)."
            )
        device = "cuda:0"
    elif device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError(
            f"Backend 'local' con device='{device}' requiere CUDA, "
            "pero no está disponible."
        )
    cfg["device"] = device

# Cliente LLM según backend
if backend == "local":
    from rdma.utils.llm_client import LocalLLMClient
    client = LocalLLMClient(
        model_type=cfg["model_type"],
        device=cfg["device"],
        cache_dir=cfg["model_cache_dir"],
    )
elif backend == "openrouter":
    from rdma.utils.llm_client import OpenRouterLLMClient
    client = OpenRouterLLMClient(model_type=cfg["model_type"])
elif backend == "api":
    from rdma.utils.llm_client import APILLMClient
    client = APILLMClient(model_type=cfg["model_type"])
elif backend == "azure":
    from rdma.utils.llm_client import AzureOpenAILLMClient
    client = AzureOpenAILLMClient(model_type=cfg["model_type"])
elif backend == "llama_cpp":
    from rdma.utils.llm_client import LlamaCppLLMClient
    client = LlamaCppLLMClient(model_type=cfg["model_type"], gguf_file=cfg.get("gguf_file"))
else:
    raise ValueError(f"Backend no soportado: {backend}")

hpo_emb = str(Path(cfg["rdma_cache_dir"]) / "vector_stores" / "G2GHPO_metadata_medembed.npy")
lab_emb = str(Path(cfg["rdma_cache_dir"]) / "tools" / "lab_tables_medembed_sm.npy")

extractor = PhenotypeExtractor(
    llm_client=client,
    extractor_type=cfg["extractor_type"],
    embeddings_file=hpo_emb,
    top_k=cfg["top_k"],
    negation=cfg["negation"],
    family_history=cfg["family_history"],
)
verifier = HPOVerifier(
    llm_client=client,
    embeddings_file=hpo_emb,
    verifier_version=cfg["verifier_version"],
    lab_embeddings_file=lab_emb,
)
matcher = HPOMatcher(
    llm_client=client,
    embeddings_file=hpo_emb,
    top_k=cfg["top_k"],
)

text = cfg["text"]
entities = extractor.extract([text])

if cfg["skip_verification"]:
    verified = entities
else:
    verified = verifier.verify(entities, text)

matched = matcher.match(verified)

# Salida JSON
result = {
    "entities": entities,
    "verified": verified if not cfg["skip_verification"] else None,
    "matched": matched,
}
print(json.dumps(result, default=str))
'''


def _write_runner(out_dir: Path) -> Path:
    """Escribe el runner como fichero temporal en output_dir."""
    runner = out_dir / "_rdma_runner.py"
    runner.write_text(_RUNNER)
    return runner


def extract_hpo(inp: ExtractHPOInput) -> ExtractHPOOutput:
    """Extrae términos HPO de un texto clínico usando RDMA.

    Ejecuta RDMA en subproceso dentro del entorno pga-phenotype,
    aislando el LLM del proceso del agente.
    """
    outdir = Path(inp.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cache = _resolve_cache_dir(inp)
    model_cache = _resolve_model_cache(inp)
    py = _python()

    runner = _write_runner(outdir)

    cfg = {
        "backend": inp.backend,
        "model_type": inp.model_type,
        "device": inp.device,
        "model_cache_dir": str(model_cache),
        "rdma_cache_dir": str(cache),
        "extractor_type": inp.extractor_type,
        "verifier_version": inp.verifier_version,
        "negation": inp.negation,
        "family_history": inp.family_history,
        "skip_verification": inp.skip_verification,
        "top_k": inp.top_k,
        "text": inp.text,
    }

    result = subprocess.run(
        [py, str(runner), json.dumps(cfg)],
        capture_output=True,
        text=True,
        cwd=str(outdir),
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"RDMA falló (exit {result.returncode}).\n"
            f"stderr:\n{result.stderr[-2000:]}"
        )

    try:
        data = json.loads(result.stdout.strip().splitlines()[-1])
    except (json.JSONDecodeError, IndexError) as e:
        raise RuntimeError(
            f"No se pudo parsear la salida de RDMA: {e}\n"
            f"stdout (últimas 20 líneas):\n"
            + "\n".join(result.stdout.splitlines()[-20:])
        )

    entities = data.get("entities", [])
    verified = data.get("verified")
    matched = data.get("matched", [])

    hpo_ids: list[str] = []
    for m in matched:
        hp = m.get("hp_id") or m.get("hpo_id")
        if hp and hp not in hpo_ids:
            hpo_ids.append(hp)

    hpo_terms_path = outdir / "hpo_terms.json"
    hpo_terms = {
        "source_text_hash": _text_hash(inp.text),
        "extractor": f"rdma:{inp.model_type}",
        "backend": inp.backend,
        "terms": [
            {
                "hpo_id": m.get("hp_id") or m.get("hpo_id"),
                "label": m.get("hpo_label") or m.get("label"),
                "matched_text": m.get("entity") or m.get("text"),
                "start": m.get("start"),
                "end": m.get("end"),
                "confidence": m.get("confidence"),
                "excluded": m.get("excluded", False),
            }
            for m in matched
        ],
        "hpo_ids": hpo_ids,
    }
    hpo_terms_path.write_text(json.dumps(hpo_terms, indent=2, default=str))

    entities_path = outdir / "entities.json"
    entities_path.write_text(json.dumps({
        "entities": entities,
        "verified": verified,
        "matched": matched,
    }, indent=2, default=str))

    counts = {
        "entities_extracted": len(entities),
        "entities_verified": len(verified) if verified else 0,
        "hpo_matched": len(matched),
        "hpo_unique": len(hpo_ids),
    }
    report_path = outdir / "extract_hpo_report.json"
    report_path.write_text(json.dumps({
        "backend": inp.backend,
        "model_type": inp.model_type,
        "extractor_type": inp.extractor_type,
        "verifier_version": inp.verifier_version,
        "skip_verification": inp.skip_verification,
        "counts": counts,
    }, indent=2))

    return ExtractHPOOutput(
        ok=True,
        tool="extract_hpo",
        message=(
            f"RDMA ({inp.backend}/{inp.model_type}): "
            f"{counts['hpo_matched']} términos HPO "
            f"(de {counts['entities_extracted']} entidades)"
        ),
        outputs={
            "hpo_terms": str(hpo_terms_path),
            "entities": str(entities_path),
            "report": str(report_path),
        },
        hpo_terms_json=hpo_terms_path,
        entities_json=entities_path,
        hpo_ids=hpo_ids,
        counts=counts,
        report_json=report_path,
    )
