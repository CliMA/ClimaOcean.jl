#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH -p sched_mit_raffaele_gpu
#SBATCH --time=120:00:00
#SBATCH --mem=200GB

module load cuda

export JULIA_CUDA_MEMORY_POOL=none

# Number of threads in SLURM mode
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}

export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

julia --project --check-bounds=no -e 'using Pkg; Pkg.update("Oceananigans")'

/nobackup1/users/ssilvest/openmpi/bin/mpiexec -np 2 julia --project --check-bounds=no prototype_omip_simulation.jl
