#!/usr/bin/env bash
# Run python from the local conda env `qchem` with Library\bin on PATH (equivalent to conda activate).
# Usage: scripts/qpy.sh script.py [args...]   or   scripts/qpy.sh -m pytest ...
ENV=/c/Users/shich/miniconda3/envs/qchem
export PATH="$ENV:$ENV/Library/bin:$ENV/Library/mingw-w64/bin:$ENV/Scripts:$PATH"
export CONDA_PREFIX="C:/Users/shich/miniconda3/envs/qchem"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export PYTHONPATH="$(cd "$(dirname "$0")/.." && pwd -W 2>/dev/null || cd "$(dirname "$0")/.." && pwd)${PYTHONPATH:+;$PYTHONPATH}"
exec "$ENV/python.exe" "$@"
