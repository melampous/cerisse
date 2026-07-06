#!/bin/bash
#SBATCH --job-name=cerisse_a100
#SBATCH -D .
#SBATCH --output=gpu_%j.out
#SBATCH --error=gpu_%j.err
#SBATCH --nodes=1                 # 8 A100 in one node (NVLink) -> best.  Use --nodes=2 if 4 GPU/node
#SBATCH --ntasks=8                # 8 MPI ranks = 8 GPUs (1 rank per GPU)
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:8              # <<< set to gpu:a100:8 if your cluster names the type
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --partition=__GPU_PARTITION__     # <<< FILL: e.g. gpu / a100 / gpuA100
#SBATCH --account=__ACCOUNT__             # <<< FILL: your allocation/project code

module purge
module load __CUDA_MODULE__ __CUDA_AWARE_MPI_MODULE__   # <<< FILL: e.g. cuda/12.2  openmpi/4.1-cuda  (or impi/2021.15)
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

# 1 rank pinned to 1 GPU; GPU-aware MPI for GPU-direct halo exchange
srun --gpu-bind=single:1 ./main3d.gnu.MPI.CUDA.ex input_bench256 amrex.use_gpu_aware_mpi=1
