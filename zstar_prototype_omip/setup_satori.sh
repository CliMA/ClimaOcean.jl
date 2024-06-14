# Upload modules
module purge all
module add spack
module add cuda/11.4
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

# MPI specific exports
export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

# Julia specific enviromental variables
export COMMON="/nobackup/users/ssilvest/perlmutter-test"
export JULIA="${COMMON}/julia/julia"

export JULIA_CUDA_MEMORY_POOL=none
export JULIA_DEPOT_PATH="${COMMON}/depot"

# Profile specific variable
export JULIA_NVTX_CALLBACKS=gc

# Number of threads in SLURM mode
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
