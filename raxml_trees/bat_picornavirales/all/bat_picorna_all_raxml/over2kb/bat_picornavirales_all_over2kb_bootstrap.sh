#!/bin/bash
#SBATCH --job-name=batpicorna_all_boot2kb
#SBATCH --partition=broadwl
#SBATCH --output=batpicorna_all_boot2kb.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa bat_picornavirales_all_over2kb_align.fasta --model TIM2+G4 --prefix T3  --seed 2 --threads auto{8}
