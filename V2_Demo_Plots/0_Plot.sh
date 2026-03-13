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

# Data preparation for plotting
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_0_Export2nc.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/V2_Demo_Plots/1_1_Data_prep.py

# ========================= Demo Fig.0
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig0a.py

# ========================= Demo Fig.1 
# Fig. 1a Exceeedance of irrigation water use, N and P runoff of the 4 major crops (wheat, rice, maize, soybean)
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig1a.py
# Fig. 1b Map of exceedance
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig1b.py
# Fig. 1c Bars of exceedance by crop
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig1c.py

# ========================= Demo Fig.2
# Fig. 2a Where should the fertilizer be reduced or increased? Map of changes in fertilizer input: (Scenario 4 - Scenario 1)/Scenario 1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig2a_data_prep.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig2a.py

# ========================= Demo Fig.3
# Fig. 3a Maps of N, P runoff reduction: (Scenario 4 - Scenario 1)/Scenario 1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig3a.py
# Fig. 3b Bar plot of changes in N, P runoff 
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig3b.py

# ========================= Demo Fig.4
# Fig. 4a Maps of crop production reduction: (Scenario 4 - Scenario 1)/Scenario 1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig4a.py
# Fig. 4b Bar plot of changes in crop production
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig4b.py

# ========================= Demo Fig.5
# Fig. 5a Plot for showing how does N runoff, N uptake, N in grain, and yield change with different fertilization reduction levels
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig5a_data_pre.py # Data preparation for sensitivity analysis
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig5a.py
# Fig. 5b P pool changes over time under different fertilization reduction levels
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig5b.py
# Fig. 5c P pool changes over time under different fertilization reduction levels
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig5c_data_pre.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig5c.py


# ========================= Statistical analyiss
# Stat1: Exceedance of irrigation water use, N and P runoff by crop and by basin
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/V2_Demo_Plots/Stat_1_Exceedance_HA.py

# Stat2: Basin avg crit N, P runoff [kg/ha] by crops
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_2_basin_avg_CritNP.py

conda deactivate