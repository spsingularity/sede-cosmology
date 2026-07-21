#!/usr/bin/env bash
# Re-run the canonical inference chain, log each driver's output.
cd "$(dirname "$0")/.."
PY=/opt/miniconda3/bin/python3
export PYTHONPATH="$PWD/src:$PWD/src/experiments"
L=results/logs
run(){ echo "=== $1 START $(date +%H:%M:%S) ==="; $PY src/experiments/$1.py "${@:2}" >$L/$1.log 2>&1; echo "=== $1 DONE rc=$? $(date +%H:%M:%S) ==="; }
run run_probe_decomposition
run run_joint_fullcmb
run run_xval_oos
run run_no_shoes_robustness
run run_tracer_loo
run run_bayesian_evidence --nlive 300
run run_barrow_mcmc --steps 800 --burn 300
run run_xval_calibration 400
echo "ALL_DONE $(date +%H:%M:%S)"
