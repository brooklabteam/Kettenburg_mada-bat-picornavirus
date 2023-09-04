#!/bin/bash
#SBATCH --job-name=batpicorna_all_1kb
#SBATCH --partition=broadwl
#SBATCH --output=batpicorna_all_1kb.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa bat_picornavirales_all_over_1kb_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{8} --bs-metric fbp,tbe
