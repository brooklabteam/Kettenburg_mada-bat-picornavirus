#!/bin/bash
#SBATCH --job-name=sapo_tree
#SBATCH --partition=broadwl
#SBATCH --output=sapo_tree.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --msa sapo_all_align.fasta --model TIM2+I+G4 --prefix T2  --seed 1 --threads auto{8}