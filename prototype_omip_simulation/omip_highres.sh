#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH -p pi_raffaele
#SBATCH --time=120:00:00
#SBATCH --mem=200GB

module load cuda

export JULIA_CUDA_MEMORY_POOL=none

# Number of threads in SLURM mode
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}

export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

/nobackup1/users/ssilvest/openmpi/bin/mpiexec -np 4 julia --project --check-bounds=no prototype_omip_simulation.jl
