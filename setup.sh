#!/usr/bin/env bash
# PublicGenomicAgent bootstrap script.
#
# Creates the base micromamba environment (pga-core), installs
# the package in editable mode, and creates the tool
# environments (pga-hts, pga-pangenome, pga-mendelian).
#
# Idempotent: safe to re-run.
#
# Requirements: micromamba in PATH.
#
# Usage:
#   ./setup.sh          # entornos completos (desarrollo local)
#   ./setup.sh --ci     # subset para CI: pga-core + pga-hts + pga-mendelian
#                       # (omite pga-pangenome, que es pesado)

set -euo pipefail

CI_MODE=0
for arg in "$@"; do
    case "$arg" in
        --ci) CI_MODE=1 ;;
        -h|--help)
            cat <<'HELP'
PublicGenomicAgent bootstrap script.

Crea los entornos micromamba del proyecto e instala el paquete
en editable mode. Idempotente.

Uso:
  ./setup.sh          # entornos completos (desarrollo local):
                      #   pga-core, pga-hts, pga-pangenome, pga-mendelian
  ./setup.sh --ci     # subset para CI:
                      #   pga-core, pga-hts, pga-mendelian
                      #   (omite pga-pangenome, que es pesado)

Requiere micromamba en PATH.
HELP
            exit 0
            ;;
        *)
            echo "Argumento desconocido: $arg" >&2
            echo "Uso: $0 [--ci]" >&2
            exit 2
            ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

if ! command -v micromamba >/dev/null 2>&1; then
    echo "ERROR: micromamba not found in PATH." >&2
    echo "Install it from https://mamba.readthedocs.io/" >&2
    exit 1
fi

PGA_ROOT="${PGA_ROOT:-$HOME/.pga}"
PGA_CORE="$PGA_ROOT/envs/pga-core"

echo "==> PGA_ROOT=$PGA_ROOT"
if [ "$CI_MODE" = "1" ]; then
    echo "==> Modo CI: solo pga-core, pga-hts, pga-mendelian"
fi
mkdir -p "$PGA_ROOT/envs"

# 1. Create pga-core if missing
if [ ! -d "$PGA_CORE/conda-meta" ]; then
    echo "==> Creating pga-core (this may take a few minutes)"
    micromamba create -y -p "$PGA_CORE" -f envs/pga-core.yml
else
    echo "==> pga-core already exists"
fi

# 2. Install the package in editable mode
echo "==> Installing publicgenomicagent in editable mode"
micromamba run -p "$PGA_CORE" pip install -e ".[dev]" -q

# 3. Create tool environments
if [ "$CI_MODE" = "1" ]; then
    TOOL_ENVS="pga-hts pga-mendelian"
else
    TOOL_ENVS="pga-hts pga-pangenome pga-mendelian"
fi

for env_name in $TOOL_ENVS; do
    if [ ! -d "$PGA_ROOT/envs/$env_name/conda-meta" ]; then
        echo "==> Bootstrapping $env_name"
        micromamba run -p "$PGA_CORE" pga env bootstrap "$env_name"
    else
        echo "==> $env_name already exists"
    fi
done

# 4. Verify
echo "==> Verifying environments"
micromamba run -p "$PGA_CORE" pga env verify || true

echo ""
echo "All set. Next steps:"
echo ""
echo "  # Local LLM (optional, requires ollama)"
echo "  ollama pull qwen2.5:3b"
echo "  $PGA_CORE/bin/pga llm ping"
echo ""
echo "  # Reproducible demo case (synthetic trio, chr2 ROI)"
echo "  $PGA_CORE/bin/python scripts/make_demo_case.py"
echo "  echo 'analiza el caso' | $PGA_CORE/bin/pga chat \\"
echo "      --case cases/demo_trio/manifest.json"
echo "  $PGA_CORE/bin/pga audit \\"
echo "      sessions/\$(ls -t sessions/ | head -n 1)/session.json --verify"
echo ""
echo "  # Local dev shortcut"
echo "  source scripts/dev-env.sh"
