#!/bin/bash
#SBATCH --job-name=pic_p2
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornavirales_refseq_P2P3_align.fasta -d nt -t ml -p 8