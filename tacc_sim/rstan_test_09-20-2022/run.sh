#!/bin/bash
#SBATCH -J sim               # Job name
#SBATCH -o sim.o             # Name of stdout output log file
#SBATCH -e sim.e             # Name of stderr output log file
#SBATCH -N 8                 # Total number of nodes to request
#SBATCH -n 64                # Total number of workers to request (distributed over nodes)
#SBATCH -p development       # The type of queue to submit to
#SBATCH -t 1:00:00           # Time limit to request (hh:mm:ss)
#SBATCH -A Extensions-to-MMREM   # Your project name
#SBATCH --mail-user=yrlee@utexas.edu # TACC will send emails with status updates
#SBATCH --mail-type=all            # Get all status updates

# load R module
module reset
module load Rstats/3.5.1

# call R code from RMPISNOW
ibrun RMPISNOW < run_sim_small_noparallel.R