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


# 1-0 Transfer the data from .csv to .nc
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_0_Export2nc.py

# 1-1 Prepare the data for Fig. 1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_1_Data_prep_Fig1.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig1.py


# # 1-2 Plot the summary plot for Fig.1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_2_Summary_Fig1_LaPlata.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_2_Summary_Fig1_Indus.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_2_Summary_Fig1_Rhine.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_2_Summary_Fig1_Yangtze.py

# Demo Plots
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_SM_Fig1_Indus.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_SM_Fig1_Yangtze.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_SM_Fig1_Rhine.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_SM_Fig1_LaPlata.py


# 2-1 Prepare the data for plotting regional boundaries of N and P runoff (sum up of all crops)
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/2_1_Data_prep_Reg_bound.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig3_1.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig3_2.py

conda deactivate