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

# ========================= Demo Fig.1 
# Fig. 1a Exceeedance of irrigation water use, N and P runoff of the 4 major crops (wheat, rice, maize, soybean)
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig1a.py
# Fig. 1b Map of exceedance
# ython /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig1b.py
# Fig. 1c Bars of exceedance by crop
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig1c.py


# ========================= Demo Fig.3
# Fig. 3a Maps of N, P runoff reduction: (Scenario 4 - Scenario 1)/Scenario 1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig3a.py
# Fig. 3b Bar plot of changes in N, P runoff 
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig3b.py


# ========================= Demo Fig.4
# Fig. 4a Maps of crop production reduction: (Scenario 4 - Scenario 1)/Scenario 1
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig4a.py
# Fig. 4b Bar plot of changes in crop production
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Demo_Fig4b.py

conda deactivate