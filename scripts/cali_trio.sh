#!/bin/bash

#SBATCH --job-name=calibration
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=50gb
#SBATCH --output=/home/yunqiyang/ukb.out
#SBATCH --error=/home/yunqiyang/ukb.err


module load gcc/12.1.0
module load R/4.2.1
# For example, Rscript calibration_sibling.R bmi "c(-1, 3)"
# Rscript calibration_trio.R bmi 
# Rscript calibration_trio.R dbp
# Rscript calibration_trio.R sbp 
Rscript calibration_trio.R diabetes 
Rscript calibration_trio.R education_yrs 

