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
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_CriticalLosses_Method1.py

# Part 1 - Appendix Fig.1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/2_Exceed.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/1_Part1_CriticalLosses.py
# -------------------------------------------------

# --------------- Plot for Part 2 -----------------
# Part 2 - Fig.2
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/3_Part2_MainPlot1_UptakeLossChange.py

# -------------------------------------------------

conda deactivate