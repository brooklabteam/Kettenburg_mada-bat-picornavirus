#!/bin/bash
#SBATCH --job-name=africa_picorna_boot3_1kb
#SBATCH --partition=broadwl
#SBATCH --output=africa_picorna_boot3_1kb.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa africa_picornaviridae_over1kb_align.fasta --model GTR+G4 --prefix T5  --seed 28 --threads auto{8}
