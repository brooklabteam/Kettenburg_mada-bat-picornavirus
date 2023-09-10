#!/bin/bash
#SBATCH --job-name=picall_ref_p2_boot
#SBATCH --partition=broadwl
#SBATCH --output=picall_ref_p2_boot.out
#SBATCH --nodes=1
#SBATCH --ntasks=
#SBATCH --ntasks-per-node=
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --bootstrap --msa .fasta --model  --prefix T3  --seed 2 --threads auto{}
