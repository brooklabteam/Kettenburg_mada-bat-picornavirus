#!/bin/bash
#SBATCH --job-name=center_tree
#SBATCH --partition=broadwl
#SBATCH --output=center_tree.out
#SBATCH --nodes=1
#SBATCH --ntasks=
#SBATCH --ntasks-per-node=
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --msa picornaviridae_refseq_all_align_center.fasta --model GTR+G4 --prefix T2  --seed 1 --threads auto{}