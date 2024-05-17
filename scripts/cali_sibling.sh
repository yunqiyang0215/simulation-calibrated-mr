#!/bin/bash

#SBATCH --job-name=calibration
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=50gb
#SBATCH --output=/home/yunqiyang/calibrated_mr/real_data_analysis/analysis/ukb.out
#SBATCH --error=/home/yunqiyang/calibrated_mr/real_data_analysis/analysis/ukb.err


module load gcc/12.1.0
module load R/4.2.1
# For example, Rscript calibration_sibling.R bmi "c(-1, 3)"
# Rscript calibration_sibling.R bmi 
# Rscript calibration_sibling.R dbp
# Rscript calibration_sibling.R sbp 

Rscript calibration_sibling.R diabetes 
Rscript calibration_sibling.R education_yrs
