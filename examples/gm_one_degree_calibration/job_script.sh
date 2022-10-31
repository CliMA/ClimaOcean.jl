#!/bin/bash

## For mpi calibration jobs we need the
## number of cores per node to be equal to
## the requested GPUs

#SBATCH -N 1
#SBATCH -n 4
#SBATCH --gres=gpu:4
#SBATCH -t 24:00:00
#SBATCH -o output.txt
#SBATCH -e error.txt
#SBATCH --mem=90GB

#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --threads-per-core=1

MPIEXECJL="/home/ssilvest/.julia/bin/mpiexecjl"

$MPIEXECJL --project -n 4 julia --check-bounds=no gm_one_degree_model_calibration.jl
