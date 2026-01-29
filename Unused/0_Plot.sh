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

# --------------- Plot for Part 1 -----------------
# Part 1 - Fig.1 
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/2_Part1_MainPlot_1_Unsus_Irri.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/2_Part1_MainPlot_2_ExceedNP.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_CriticalLosses_Method1.py

# Part 1 - Appendix Fig.1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/2_Exceed.py


# Part 2 - Plot excessive N, P losses
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_Current_Losses_Total.py


# ---------------- Testing plot -----------
# Unit [kg/ha]
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_CriticalLosses_Method1.py
# Unit [kg]
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_CriticalLosses_Method1_Total.py
# Unit [kg/ha]
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_CriticalLosses_Method2.py
# Unit [kg]
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_CriticalLosses_Method2_Total.py
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Unused/1_Part1_CriticalLosses_Method3.py
# -------------------------------------------------

# --------------- Plot for Part 2 -----------------
# Part 2 - Fig.2
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/3_Part2_MainPlot1_UptakeLossChange.py
# -------------------------------------------------


# --------------- Plot for Part 3 -----------------
# Part 3 - Main figure 3
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/4_Part3_MainPlot3_ProDecrea.py
# -------------------------------------------------



conda deactivate