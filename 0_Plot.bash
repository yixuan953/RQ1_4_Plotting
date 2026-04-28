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

# ================== Data preparation
# 1. Trasnform the .csv format to .nc

# 2. Summarize the simulated output 
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Data_Prep1_csv2nc.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Data_Prep2_summary.py
# ----> Results were saved in "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary"

# # 3. Check where the boundaries have been exceeded: To plot Fig 3b
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig3b_BoundaryCheck.py

# # ================== Statistical analysis
# # 1-1. Calculate the boundaries for total N and P delivery, agricultural runoff, cropland runoff, and each crop types [ktons]
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig2a_1.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig2a_2.py

# # 1-2. Calculate the average crop-specific N and P runoff to surface water [kg/ha]
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Table1_Basin_avg_CritNP.py
# # ---> Results were saved in "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Statistics/1_Boundary"

# # 2. Calculate the simulated crop production and runoff
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig3a_1.py

# # 3. Calculate the share of harvested area where the boundaries have been exceeded
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig4c_Exceedance_HA.py

# # 4. Calculate how much fertilizers should be reduced
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig5a_Fert_Red.py

# # 5. Calculate the crop yield and N, P runoff reduction after fertilizer and (or) irrigation reduction
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig5_Implications.py

# ================== Plotting
# 1. Plot the barplots and maps for the quantified regional safe boundaries
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig2a_BarPlorts.py
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig2b_Maps.py

# 2. Plot the simulated crop production and runoff
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig3a_SimRunoffProd_bars.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig3b_SimRunoffProd_maps.py

# 3. Compare the simulated runoff with the boundaries
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig4a_BarPlots.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig4b_BoundaryCheck.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig4c_ExceedanceMaps.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig4c_Exceedance_bars.py

# 5. Plot the crop yield and N, P runoff reduction after fertilizer reduction
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig5a_Fert_Rate_Red.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig5b_ExceedanceMaps.py
# python  /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig5c_Yield_Red_Map.py

conda deactivate