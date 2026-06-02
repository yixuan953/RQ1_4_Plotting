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

# ================== Data preparation =====================
# 1. Trasnform the .csv format to .nc
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Data_Prep1_csv2nc.py

# 2. Summarize the simulated output 
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Data_Prep2_summary.py

# ===================== Main figure statsitics and plotting =====================
# --------->  Main Figure 2: Boundary exceedance under the current irrigation and fertilization practices

# Fig 2a & 2b - Bar plot of actual vs critical runoff for each crop type and basin
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig2a_2.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig2ab_BarPlots_v2.py


# Fig 2c - Areas with boundaries' exceedances & Share of harvested area exceeding sustainable irrigation amount or N (P) runoff boundaries
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig2c_BoundaryCheck.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig2c_2.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig2c_Exceedance_bars.py

# Fig 2d & 2e - Maps of boundaries' exceedance [kg N/ha/yr or kg P/ha/yr] for each crop type and basin
python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig2d_ExceedanceMaps.py

# --------->  Main Figure 3: Impacts on crop production when reducing irrigation and fertilization
# Fig 3a - Crop production under different scenarios
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig3a.py
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig3a_CropProd_bars.py

# Fig 3b - Crop production reduction [%] with sustainable irrigation and feritlization
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig3b_CropProd_red_Maps.py

# --------->  Main Figure 4: Nutrient losses under different scenarios
# Fig 4a - Nutrient losses under different scenarios 
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig4a_NP_runoff_bars.py

# Fig 4b - Maps of N and P runoff exceedance with sustainable irrigation and reduced fertilization
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig4b_NP_exceedance_Maps.py

# ---------> Main Figure 5: Trade-offs between crop production and nutrient losses under different scenarios
# Fig 5a ---> Basin-scale crop production vs. N, P runoff when allowing different percentage of yield reduction
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig5a_1.py # Basin-scale crop production vs. N, P runoff when allowing different percentage of yield reduction
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig5a_2.py # Boundaries without contribution from sewage
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Stat_Fig5a_3.py # Boundaries without contribution from sewage
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Fig5a_Prod_Losses_tradeoffs.py


# ---------> Extended Figures:
# Extended Fig. 1 Crop-specific regional safe boundaries for N runoff, simulated current crop-specific N runoff, and the boundaries’ exceedance.
# Extended Fig. 2 Crop-specific regional safe boundaries for P runoff, simulated current crop-specific P runoff, and the boundaries’ exceedance.
# Extended Fig. 3 Simulated major crops’ production with current (2015) irrigation and fertilization
# Extended Fig. 4 Simulated crop-specific N runoff after reducing irrigation and fertilization and the boundaries’ exceedance.
# Extended Fig. 5 Simulated crop-specific P runoff after reducing irrigation and fertilization and the boundaries’ exceedance.
# Extended Fig. 6 Simulated major crops’ production with reduced irrigation and fertilization
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Extended_Fig5.py


# ---------> Sensitivity analysis:
# python /lustre/nobackup/WUR/ESG/zhou111/1_RQ1_Code/4_Plotting/Sens_Fig1.py

conda deactivate