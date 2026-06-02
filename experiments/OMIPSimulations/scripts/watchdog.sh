#!/bin/bash
# Watchdog that keeps store.sh jobs alive for the given run names.
# Usage: ./watchdog.sh orca_ncar orca_corrected_snow_cb0.1
# Run inside tmux from the same directory as store.sh.

set -euo pipefail

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <run_name1> [run_name2] ..."
    echo "Example: $0 orca_ncar orca_corrected_snow_cb0.1"
    exit 1
fi

RUN_NAMES=("$@")

while true; do
    for run in "${RUN_NAMES[@]}"; do
        if ! squeue -u "$USER" -n "store_${run}" -h | grep -q .; then
            echo "$(date): store_${run} not found, relaunching"
            ./store.sh "$run"
        else
            echo "$(date): store_${run} is running"
        fi
    done
    sleep 3600
done
