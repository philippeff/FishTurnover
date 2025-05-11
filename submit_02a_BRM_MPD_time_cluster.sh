#!/bin/bash
#SBATCH --mail-user=<pfa15@sfu.ca>
#SBATCH --mail-type=ALL
#SBATCH --job-name=BRM_MPD_time
#SBATCH --account=def-mooers
#SBATCH --output=./out/%x-%j.out
#SBATCH --cpus-per-task=4		# number of processes ("ncores" for parallelization in R)
#SBATCH --nodes=1             	# number of node MUST be 1
#SBATCH --mem-per-cpu=8Gb    	# memory; default unit is megabytes
#SBATCH --time=7-00:00:00      	# time (DD-HH:MM:SS). Béluga, Graham and Narval = max 7 days. Cedar = max 28 days.

#set the environment
module load StdEnv/2020
module load r/4.2.2

Rscript 02a_BRM_MPD_time_cluster.R