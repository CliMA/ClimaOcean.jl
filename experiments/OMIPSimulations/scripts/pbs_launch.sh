#!/bin/bash

#PBS -A UMIT0080
#PBS -N tenthdegree
#PBS -j oe
#PBS -q main
#PBS -l walltime=12:00:00
#PBS -o tenthdegree.out
#PBS -l select=2:ncpus=64:mpiprocs=4:ngpus=4:gpu_type=a100:mem=384GB

# Submit a tenth-degree OMIP simulation on Derecho (PBS + cray-mpich).
#
# Run from this directory (`experiments/OMIPSimulations/scripts/`):
#     qsub pbs_launch.sh
#
# Override defaults with env vars at submit time, e.g.:
#     qsub -v RUN_NAME=tenthdegree_corrected pbs_launch.sh
#
# For halfdegree/ORCA configurations use `launch.sh` (SLURM).

set -euo pipefail

export JULIA_DEPOT_PATH=/glade/work/$USER/.julia
export JULIA_NUM_PRECOMPILE_TASKS=64
export JULIA_NUM_THREADS=64

module --force purge
module load ncarenv nvhpc cuda cray-mpich

export MPICH_GPU_SUPPORT_ENABLED=1
export JULIA_MPI_HAS_CUDA=true
export PALS_TRANSFER=false
export JULIA_CUDA_MEMORY_POOL=none

# ── Inputs ───────────────────────────────────────────────────────────────
RUN_NAME="${RUN_NAME:-tenthdegree}"
JULIA="${JULIA:-$HOME/julia-1.12.5/bin/julia}"
JULIA_THREADS="${JULIA_THREADS:-16}"

# ── Julia simulation expression (tenth-degree on 1×8 distributed GPUs) ────
# `forcing_dir` is set to the JRA55 default scratch cache, so the JRA55
# data is read in place — no separate staging directory.
JULIA_EXPR="using OMIPSimulations
using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations
using NumericalEarth
using CUDA

sim = omip_simulation(:tenthdegree;
                      arch = Distributed(GPU(), partition=Partition(1, 8)),
                      Nz = 100,
                      depth = 5500,
                      Δz_top = 2.5,
                      κ_skew = nothing,
                      κ_symmetric = nothing,
                      biharmonic_timescale = nothing,
                      Δt = 2minutes,
                      forcing_dir = NumericalEarth.DataWrangling.JRA55.download_JRA55_cache,
                      file_splitting_interval = 180days,
                      output_dir = \"${RUN_NAME}_run\",
                      filename_prefix = \"${RUN_NAME}\")

sim.stop_time = 91days
run!(sim)

sim.Δt = 10minutes
sim.stop_time = 300 * 365days
run!(sim; pickup = true)"

mpiexec -n 8 -ppn 4 "$JULIA" --project=.. --check-bounds=no -t "${JULIA_THREADS}" -e "$JULIA_EXPR"
