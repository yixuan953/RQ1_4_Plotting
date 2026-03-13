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


# ========================= Demo Fig.1: Boundaries for cropland N and P runoff to surface water in ktons & kg/ha
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig1_Stat1.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig1_Stat2.py
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig1a_BoundariesBar.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig1b_BoundariesMaps.py

# ========================= Demo Fig.2: Simulated yield and losses
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig2_Yield_Losses.py

# ========================= Demo Fig.3: Exceedance of boundaries
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig3_Exceedance.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig3_BarPlots.py


# ========================= Statistical analyiss


conda deactivate