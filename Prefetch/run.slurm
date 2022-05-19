#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH -o job.%N.%j.out  # STDOUT
#SBATCH -e job.%N.%j.err  # STDERR
#SBATCH --job-name prefetch
#SBATCH --exclusive

module purge
module load userspace/custom opt/all userspace/all
module load cmake/3.22.2
module load intel-compiler/64/2018.3.222
module load intel-mkl/64/2018.3.222
module load intel-mpi/64/2018.3.222
module load intel-mpi/64/2018.3.222

make all
