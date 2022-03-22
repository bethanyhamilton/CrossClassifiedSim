#!/bin/bash
#SBATCH -J CC                  # Job name
#SBATCH -o CC.o%j              # Name of stdout output file (%j expands to jobId)
#SBATCH -e CC.o%j              # Name of stderr output file(%j expands to jobId)
#SBATCH -p development          # Submit to the 'normal' or 'development' queue
#SBATCH -N 16                    # Total number of nodes ??
#SBATCH -n 272                  # Total number of mpi tasks requested ???
#SBATCH -t 2:00:00             # Run time (hh:mm:ss)
#SBATCH --mail-user=bethanyhamilton@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH -Cross Classified

# load R module -- change
module load Rstats/3.5.1

# call R code from RMPISNOW
ibrun RMPISNOW < ./run_sim.R
