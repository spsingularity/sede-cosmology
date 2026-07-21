#!/usr/bin/env bash
cd "$(dirname "$0")/.."
export PYTHONPATH="$PWD/src:$PWD/src/experiments"
export PYTHONUNBUFFERED=1
echo "START $(date)"
# caffeinate -i keeps the system awake (no idle sleep) for the duration of this run
caffeinate -i /opt/miniconda3/bin/python3 -u src/experiments/run_xval_calibration.py 2000 200
echo "END rc=$? $(date)"
