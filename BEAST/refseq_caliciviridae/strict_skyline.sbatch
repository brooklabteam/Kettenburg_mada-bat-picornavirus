#!/bin/sh
#SBATCH --job-name=strict_beast_cal
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=8	#number of tasks
#SBATCH --cpus-per-task=1	#number of threads per task
#SBATCH --time=36:00:00

module load beagle/trunk
module load beast/2.6.2

/software/beast-2.6.2-el7-x86_64/bin/beast -beagle_cpu -seed 777 strict_skyline.xml 

#/software/beast-2.6.2-el7-x86_64/bin/beast -beagle_cpu -seed 88 -resume strict_skyline.xml 
