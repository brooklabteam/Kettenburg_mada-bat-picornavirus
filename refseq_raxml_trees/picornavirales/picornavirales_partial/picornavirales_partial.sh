#!/bin/bash
#SBATCH --job-name=picornavirales_partial
#SBATCH --partition=broadwl
#SBATCH --output=picornavirales_partial.out
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornavirales_partial_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{12} --bs-metric fbp,tbe