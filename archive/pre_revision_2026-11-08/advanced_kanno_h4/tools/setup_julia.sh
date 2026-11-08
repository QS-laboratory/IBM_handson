#!/usr/bin/env bash
set -euo pipefail
H4_EXERCISE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export JULIA_DEPOT_PATH="$H4_EXERCISE_ROOT/.julia_depot"
export PYTHON="$(command -v python3)"
export PYTHONPATH="$H4_EXERCISE_ROOT/.python_deps${PYTHONPATH:+:$PYTHONPATH}"
export JULIA_NUM_THREADS=2
julia --startup-file=no --project="$H4_EXERCISE_ROOT/julia" -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall"); Pkg.status()'
