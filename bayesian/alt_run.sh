#!/bin/bash
#SBATCH --time=0-08:00:00
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=outfile
module purge
module load R

### Run simulation code and benchmark the start/end times
start_time=$(date +%s)
Rscript mcmc_divided.R
finish_time=$(date +%s)

### Calculate and output total runtime
elapsed_time=$((finish_time  - start_time))
((sec=elapsed_time%60, elapsed_time/=60, min=elapsed_time%60, hrs=elapsed_time/60))
timestamp=$(printf "Total time taken - %d hours, %d minutes, and %d seconds." $hrs $min $sec)
echo $timestamp

