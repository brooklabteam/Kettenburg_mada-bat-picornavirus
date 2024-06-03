#!/bin/bash
#SBATCH --job-name=center_boot2
#SBATCH --partition=broadwl
#SBATCH --output=center_boot2.out
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa picornaviridae_refseq_all_align_center.fasta --model GTR+G4 --prefix T4  --seed 2 --threads auto{10}
