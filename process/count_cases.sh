#!/bin/bash
#
#SBATCH --job-name=count_cases
#SBATCH --output=log_count_cases.txt

#SBATCH --mem=64G


module load R/3.5.0

Rscript ./counts_cases.R
