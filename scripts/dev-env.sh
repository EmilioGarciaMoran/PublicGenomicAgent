#!/usr/bin/env bash
export MAMBA_ROOT_PREFIX=~/.pga
eval "$(micromamba shell hook -s bash)"
micromamba activate pga-core
