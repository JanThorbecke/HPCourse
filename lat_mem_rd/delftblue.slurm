#!/bin/bash
#
#SBATCH --job-name="lat_mem_rd"
#SBATCH --partition=compute
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=20G

module use /mnt/shared/apps/modulefiles
module load 2022r1
module load compute

srun ./lat_mem_rd 16384 128 >& delftblue

