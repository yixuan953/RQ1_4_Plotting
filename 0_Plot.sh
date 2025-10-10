#!/bin/bash
#-----------------------------Mail address-----------------------------

#-----------------------------Output files-----------------------------
#SBATCH --output=HPCReport/output_%j.txt
#SBATCH --error=HPCReport/error_output_%j.txt

#-----------------------------Required resources-----------------------
#SBATCH --time=600
#SBATCH --mem=250000

#--------------------Environment, Operations and Job steps-------------
source /home/WUR/zhou111/miniconda3/etc/profile.d/conda.sh

conda activate myenv

# Part 1: Get the yield and exceedance of N, P losses
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_Yield_Losses.py 

conda deactivate