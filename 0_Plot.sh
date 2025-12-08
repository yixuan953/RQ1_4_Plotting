#!/bin/bash
#-----------------------------Mail address-----------------------------

#-----------------------------Output files-----------------------------
#SBATCH --output=HPCReport/output_%j.txt
#SBATCH --error=HPCReport/error_output_%j.txt

#-----------------------------Required resources-----------------------
#SBATCH --time=60
#SBATCH --mem=25000

#--------------------Environment, Operations and Job steps-------------
source /home/WUR/zhou111/miniconda3/etc/profile.d/conda.sh
conda activate myenv


# 1-1 Prepare the data for Fig. 1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_1_Data_prep_Fig1.py

# 1-2 Plot the summary plot for Fig.1
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_2_Summary_Fig1.py

conda deactivate