#!/bin/bash
#SBATCH --job-name=batpic_all_p3_boot
#SBATCH --partition=broadwl
#SBATCH --output=batpic_all_p3_boot.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa bat_picornavirales_all_P1_align.fasta --model GTR+G4 --prefix T5  --seed 9 --threads auto{8}