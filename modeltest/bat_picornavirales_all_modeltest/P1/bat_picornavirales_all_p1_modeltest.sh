#!/bin/bash
#SBATCH --job-name=bat_pic_all_p1
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i bat_picornavirales_all_P1_align.fasta -d nt -t ml -p 8
