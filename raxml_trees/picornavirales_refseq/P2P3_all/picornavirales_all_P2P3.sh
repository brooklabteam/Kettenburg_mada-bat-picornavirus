#!/bin/bash
#SBATCH --job-name=picall_ref_p2
#SBATCH --partition=broadwl
#SBATCH --output=picall_ref_p2.out
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornavirales_all_P2.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{12} --bs-metric fbp,tbe
