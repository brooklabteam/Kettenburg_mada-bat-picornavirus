#!/bin/bash
#SBATCH --job-name=picornaviridae_all_boot4
#SBATCH --partition=broadwl
#SBATCH --output=picornaviridae_all_boot4.out
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa picornaviridae_all_align.fasta --model TVM+G4 --prefix T6  --seed 5 --threads auto{10}
