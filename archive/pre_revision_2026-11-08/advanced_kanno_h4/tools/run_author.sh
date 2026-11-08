#!/usr/bin/env bash
set -euo pipefail
H4_EXERCISE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export JULIA_DEPOT_PATH="$H4_EXERCISE_ROOT/.julia_depot"
export PYTHONPATH="$H4_EXERCISE_ROOT/.python_deps${PYTHONPATH:+:$PYTHONPATH}"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
export PYTHON="$(command -v python3)"
export H4_GEOMETRY="${H4_GEOMETRY:-0}"
export H4_OPT_SWEEPS="${H4_OPT_SWEEPS:-1000}"
if [[ $# -eq 0 ]]; then set -- canonical boys_ordered boys; fi
julia --startup-file=no --project="$H4_EXERCISE_ROOT/julia" "$H4_EXERCISE_ROOT/julia/run_author.jl" "$@"
