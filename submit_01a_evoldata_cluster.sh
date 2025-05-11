#!/bin/bash
#SBATCH --mail-user=<pfa15@sfu.ca>
#SBATCH --mail-type=ALL
#SBATCH --job-name=evoldata
#SBATCH --account=def-mooers
#SBATCH --output=./out/%x-%j.out
#SBATCH --nodes=1             	# number of node MUST be 1
#SBATCH --cpus-per-task=32		# number of processes ("ncores" for parallelization in R)
#SBATCH --mem-per-cpu=4Gb    	# memory; default unit is megabytes
#SBATCH --time=1-00:00:00      # time (DD-HH:MM:SS)

#set the environment
module load StdEnv/2020
module load r/4.2.2

Rscript 01a_evoldata_cluster.R