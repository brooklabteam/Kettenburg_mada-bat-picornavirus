#!/bin/bash
#SBATCH --job-name=picornaviridae_partial_boot
#SBATCH --partition=broadwl
#SBATCH --output=picornaviridae_partial_boot.out
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa picornaviridae_partial_align.fasta --model GTR+G4 --prefix T3  --seed 1 --threads auto{10}
