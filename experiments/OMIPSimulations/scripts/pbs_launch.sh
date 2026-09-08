#!/bin/bash

#PBS -A UMIT0080
#PBS -N twelfthdegree
#PBS -j oe
#PBS -q main
#PBS -l walltime=12:00:00
#PBS -o twelfthdegree.out
#PBS -l select=2:ncpus=64:mpiprocs=4:ngpus=4:gpu_type=a100:mem=384GB

# Submit a tenth-degree OMIP simulation on Derecho (PBS + cray-mpich).
#
# Run from this directory (`experiments/OMIPSimulations/scripts/`):
#     qsub pbs_launch.sh
#
# Override defaults with env vars at submit time, e.g.:
#     qsub -v RUN_NAME=twelfthdegree_corrected pbs_launch.sh
#
# For halfdegree/ORCA configurations use `launch.sh` (SLURM).

set -euo pipefail

export JULIA_DEPOT_PATH=/glade/work/$USER/.julia
export JULIA_NUM_PRECOMPILE_TASKS=64
export JULIA_NUM_THREADS=64

module --force purge
# cray-mpich is hierarchical on Derecho: it needs ncarenv + a compiler (nvhpc) loaded first.
# Pin exact versions so the batch stack matches the one CUDA.jl was set_runtime_version!'d against
# (local CUDA 12.9 toolkit; see LocalPreferences.toml in the project).
module load ncarenv/25.10 nvhpc/26.1 cuda/12.9.0 cray-mpich/8.1.32

export MPICH_GPU_SUPPORT_ENABLED=1
export JULIA_MPI_HAS_CUDA=true
export PALS_TRANSFER=false
export JULIA_CUDA_MEMORY_POOL=none

# cray-mpich dlopen'd by MPI.jl does not auto-link the GPU Transport Layer, so with
# MPICH_GPU_SUPPORT_ENABLED=1 it aborts ("GTL library is not linked") unless we preload it.
export LD_PRELOAD="$CRAY_MPICH_ROOTDIR/gtl/lib/libmpi_gtl_cuda.so${LD_PRELOAD:+:$LD_PRELOAD}"

# The PBS prologue sets CUDA_VISIBLE_DEVICES to the MOM node's GPU *UUIDs*, and mpiexec forwards
# that single value to every rank on every node. On any node other than the MOM, those UUIDs
# do not exist, so CUDA inits with NO_DEVICE and GPU() throws "a CUDA GPU was not found". Unset it
# so each node sees its own 4 local GPUs (cgroup-restricted); Oceananigans then binds each
# node-local rank to a distinct device via the COMM_TYPE_SHARED communicator.
unset CUDA_VISIBLE_DEVICES

# ── Inputs ───────────────────────────────────────────────────────────────
RUN_NAME="${RUN_NAME:-twelfthdegree}"
JULIA="${JULIA:-$HOME/software/julia-1.12.6/bin/julia}"
JULIA_THREADS="${JULIA_THREADS:-16}"

# ── Julia simulation expression (tenth-degree on 1×8 distributed GPUs) ────
# `forcing_dir` is set to the JRA55 default scratch cache, so the JRA55
# data is read in place — no separate staging directory.
JULIA_EXPR="using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations
using NumericalEarth
using CUDA

sim = omip_simulation(:twelfthdegree;
                      arch = Distributed(GPU(), partition=Partition(1, 8)),
                      Nz = 100,
                      depth = 5500,
                      Δz_top = 1.5,
                      κ_skew = nothing,
                      κ_symmetric = nothing,
                      biharmonic_timescale = nothing,
                      Δt = 5minutes,
                      forcing_dir = NumericalEarth.DataWrangling.JRA55.download_JRA55_cache,
                      checkpoint_interval = 180days,
                      file_splitting_interval = 180days,
                      output_dir = \"${RUN_NAME}_run\",
                      filename_prefix = \"${RUN_NAME}\")

# Spin up at the safe 5-minute step until just past the first checkpoint (180 days),
# so a checkpoint is on disk before switching to the 10-minute step. This also lets
# resubmissions resume via pickup instead of redoing the spin-up from scratch.
sim.stop_time = 181days
run!(sim)

sim.Δt = 10minutes
sim.stop_time = 300 * 365days
# Pick up the furthest-progressed checkpoint by iteration, NOT by mtime: a resubmission
# re-runs the spin-up and rewrites the day-90 checkpoint (cleanup=false keeps all files),
# so `pickup=true` (:recent_time_stamp) would rewind to day 90. :highest_iteration avoids that.
run!(sim; pickup = :highest_iteration)"

mpiexec -n 8 -ppn 4 "$JULIA" --project=.. --check-bounds=no -t "${JULIA_THREADS}" -e "$JULIA_EXPR"
