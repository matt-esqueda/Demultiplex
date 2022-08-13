#!/bin/bash

#SBATCH --partition=bgmp                    ### Partition (like a queue in PBS)
#SBATCH --job-name=demux                    ### Job Name
#SBATCH --nodes=1                           ### Number of nodes needed for the job
#SBATCH --cpus-per-task=1                   ### Number of cpus per task
#SBATCH --account=bgmp                      ### Account used for job submission
#SBATCH --output=reads-%j.out               ### File to store output
#SBATCH --error=reads-%j.error              ### File to store error

/usr/bin/time -v

python ./demultiplex.py -ic 30