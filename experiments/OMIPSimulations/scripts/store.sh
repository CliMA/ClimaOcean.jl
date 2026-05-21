#!/bin/bash
# Move completed OMIP outputs from one or more live run folders to
# $DATA/OMIP-data/<RUN_NAME>_run while the launch.sh jobs are still
# running. Multiple run names can be passed; the same sbatch job cycles
# through all of them once per hour.
#
# Logic (per run):
#   - Part files (*_part<N>.jld2): the highest N per filename group is
#     still being written by the running sim, so it is left in place;
#     all older parts are moved.
#   - Checkpoint files (*_checkpoint_iteration<N>.jld2): the highest
#     iteration is kept locally so `run!(sim; pickup=true)` still works;
#     older checkpoints are moved.
#   - Anything else in the run folder is left untouched.
#   - If a run directory disappears mid-cycle (sim crashed, folder
#     deleted), it is skipped with a warning rather than killing the job.
#
# Must be run from the same directory as launch.sh (i.e. this scripts
# folder) so that <RUN_NAME>_run resolves the same way it does for the
# running simulations.
#
# Usage:
#   ./store.sh orca
#   ./store.sh orca_ncar
#   ./store.sh orca_corrected_snow_cb0.1
#   ./store.sh orca_nori orca_rbvd orca_catke
#
# The arguments are RUN_NAMEs (same as the job names from launch.sh).
# DATA must be set in the calling shell (it is propagated to the
# sbatch job via --export=ALL).

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: ./store.sh <run_name> [<run_name> ...]

The <run_name>s match the RUN_NAME built by launch.sh, e.g.:
  orca, orca_ncar, orca_corrected_snow, orca_ncar_cb0.1, halfdegree, ...

Examples:
  ./store.sh orca
  ./store.sh orca_ncar
  ./store.sh orca_corrected_snow_cb0.1
  ./store.sh orca_nori orca_rbvd orca_catke
USAGE
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

case "${1:-}" in
    -h|--help)
        usage
        exit 0
        ;;
esac

RUN_NAMES=("$@")

if [[ -z "${DATA:-}" ]]; then
    echo "Error: DATA environment variable is not set" >&2
    exit 1
fi

# Validate all run directories exist up-front so submission fails fast.
for RUN_NAME in "${RUN_NAMES[@]}"; do
    RUN_DIR="${RUN_NAME}_run"
    if [[ ! -d "$RUN_DIR" ]]; then
        echo "Error: run directory '$RUN_DIR' not found in $(pwd)" >&2
        echo "       (store.sh must be run from the same directory as launch.sh)" >&2
        exit 1
    fi
done

# Tag for output files / job name: single run keeps the old naming,
# multi-run collapses to a count-tagged identifier.
if [[ ${#RUN_NAMES[@]} -eq 1 ]]; then
    TAG="${RUN_NAMES[0]}"
else
    TAG="multi${#RUN_NAMES[@]}"
fi

JOB_NAME="${JOB_NAME:-store_${TAG}}"

SBATCH_ARGS=()
SBATCH_ARGS+=(-o "store_${TAG}.out")
SBATCH_ARGS+=(-e "store_${TAG}.err")
SBATCH_ARGS+=(-J "$JOB_NAME")
SBATCH_ARGS+=(--export="ALL,RUN_NAMES_STR=${RUN_NAMES[*]}")

sbatch "${SBATCH_ARGS[@]}" <<'EOF'
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p sched_mit_raffaele
#SBATCH --time=24:00:00
#SBATCH --mem=4GB

set -uo pipefail

read -ra RUN_NAMES <<< "$RUN_NAMES_STR"

echo "Storing outputs for ${#RUN_NAMES[@]} run(s): ${RUN_NAMES[*]}"

shopt -s nullglob

archive_run() {
    local RUN_NAME="$1"
    local RUN_DIR="${RUN_NAME}_run"
    local DEST_DIR="${DATA}/OMIP-data/${RUN_DIR}"

    echo ""
    echo "=== ${RUN_NAME} ==="
    echo "  source: $(pwd)/${RUN_DIR}"
    echo "  dest:   ${DEST_DIR}"

    if [[ ! -d "$RUN_DIR" ]]; then
        echo "  warning: run directory '$RUN_DIR' missing — skipping this cycle" >&2
        return 0
    fi

    mkdir -p "$DEST_DIR"

    # ------------------------------------------------------------------
    # Part files: *_part<N>.jld2
    # The highest N per filename group is still being written, so it is
    # left in place; everything older is moved.
    # ------------------------------------------------------------------
    local -A max_part
    local f base tail n group current max
    for f in "$RUN_DIR"/*_part[0-9]*.jld2; do
        base=$(basename "$f")
        tail="${base##*_part}"
        n="${tail%.jld2}"
        [[ "$n" =~ ^[0-9]+$ ]] || continue
        group="${base%_part${n}.jld2}"
        current="${max_part[$group]:-0}"
        if (( n > current )); then
            max_part[$group]=$n
        fi
    done

    local moved_parts=0 kept_parts=0
    for f in "$RUN_DIR"/*_part[0-9]*.jld2; do
        base=$(basename "$f")
        tail="${base##*_part}"
        n="${tail%.jld2}"
        [[ "$n" =~ ^[0-9]+$ ]] || continue
        group="${base%_part${n}.jld2}"
        max="${max_part[$group]:-0}"
        if (( n == max )); then
            echo "  skip (active):  ${base}"
            kept_parts=$((kept_parts + 1))
            continue
        fi
        echo "  move:           ${base}"
        mv -- "$f" "$DEST_DIR/"
        moved_parts=$((moved_parts + 1))
    done

    # ------------------------------------------------------------------
    # Checkpoint files: *_iteration<N>.jld2
    # The latest iteration per group is required for run!(sim; pickup=true)
    # so it is kept locally; earlier checkpoints are moved.
    # ------------------------------------------------------------------
    local -A max_ckpt
    for f in "$RUN_DIR"/*_iteration[0-9]*.jld2; do
        base=$(basename "$f")
        tail="${base##*_iteration}"
        n="${tail%.jld2}"
        [[ "$n" =~ ^[0-9]+$ ]] || continue
        group="${base%_iteration${n}.jld2}"
        current="${max_ckpt[$group]:-0}"
        if (( n > current )); then
            max_ckpt[$group]=$n
        fi
    done

    local moved_ckpts=0 kept_ckpts=0
    for f in "$RUN_DIR"/*_iteration[0-9]*.jld2; do
        base=$(basename "$f")
        tail="${base##*_iteration}"
        n="${tail%.jld2}"
        [[ "$n" =~ ^[0-9]+$ ]] || continue
        group="${base%_iteration${n}.jld2}"
        max="${max_ckpt[$group]:-0}"
        if (( n == max )); then
            echo "  skip (latest):  ${base}"
            kept_ckpts=$((kept_ckpts + 1))
            continue
        fi
        echo "  move:           ${base}"
        mv -- "$f" "$DEST_DIR/"
        moved_ckpts=$((moved_ckpts + 1))
    done

    echo "  done: moved ${moved_parts} part file(s) (kept ${kept_parts})," \
         "moved ${moved_ckpts} checkpoint file(s) (kept ${kept_ckpts})."
}

while true; do
    for RUN_NAME in "${RUN_NAMES[@]}"; do
        archive_run "$RUN_NAME"
    done

    echo ""
    echo "Sleeping for 1 hour"
    sleep 3600
done
EOF
