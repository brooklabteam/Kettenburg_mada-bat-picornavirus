#!/bin/bash
#SBATCH --job-name=cal_mod
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i calici_refseq_all_align.fasta -d nt -t ml -p 8
