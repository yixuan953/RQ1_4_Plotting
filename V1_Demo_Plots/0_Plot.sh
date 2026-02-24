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
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/V1_Demo_Plots/1_0_Export2nc.py

# 1-1 Prepare the data for Fig. 1
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/V1_Demo_Plots/1_1_Data_prep_Fig1.py
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
# 2-2 Prepare the data for plotting fertilizer reducton of N and P runoff (sum up of all crops)
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/2_2_Data_prep_Fert_red.py

# Demo Plots: Crop production
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig2_1.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig2_2.py

# Demo Plots: N, P exceedance
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig3_1.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig3_2.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_SM_Fig3.py

# Demo Plots: annual N, P losses changes
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig4_1.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig4_2_1.py # Runoff
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig4_2_2.py # P pool
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig4_3.py

# Demo Plots: NP interactions
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Main_Fig5.py

conda deactivate