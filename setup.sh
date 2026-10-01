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

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if ! command -v micromamba >/dev/null 2>&1; then
    echo "ERROR: micromamba not found in PATH." >&2
    echo "Install it from https://mamba.readthedocs.io/" >&2
    exit 1
fi

PGA_ROOT="${PGA_ROOT:-$HOME/.pga}"
PGA_CORE="$PGA_ROOT/envs/pga-core"

echo "==> PGA_ROOT=$PGA_ROOT"
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
for env_name in pga-hts pga-pangenome pga-mendelian; do
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
echo "  # End-to-end demo on a synthetic recessive trio"
echo "  python3 scripts/make_demo_trio.py"
echo "  $PGA_CORE/bin/pga run-trio \\"
echo "      --father results/demo_trio/father.bam \\"
echo "      --mother results/demo_trio/mother.bam \\"
echo "      --proband results/demo_trio/proband.bam \\"
echo "      --reference results/demo_trio/ref.fa \\"
echo "      --region chr1:1000-1200 \\"
echo "      --case-id DEMO_TRIO --label DEMO1"

