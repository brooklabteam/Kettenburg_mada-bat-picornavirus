#!/bin/bash
#SBATCH --job-name=full1kb
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i bat_picornavirales_full_over1kb_align.fasta -d nt -t ml -p 8
