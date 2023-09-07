#!/bin/bash
#SBATCH --job-name=ifla_boot2
#SBATCH --partition=broadwl
#SBATCH --output=ifla_boot2.out
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --ntasks-per-node=9
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa iflaviridae_refseq_align.fasta --model GTR+G4 --prefix T4  --seed 21 --threads auto{9}
