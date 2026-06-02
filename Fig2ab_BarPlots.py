import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker

# --- Configuration (Maintained) ---
Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
InputCrops = ["soybean", "mainrice", "secondrice", "maize", "winterwheat"]
FinalCategories = ["Soybean", "Rice", "Maize", "Wheat"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig2"

# input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
# output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig_Extra"

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"

crop_colors = {'Wheat': "#EBBB7C", 'Maize': "#B36767", 'Rice': "#C89CDC", 'Soybean': "#42949F"}
crop_mapper = {"soybean": "Soybean",  "mainrice": "Rice", "secondrice": "Rice","maize": "Maize", "winterwheat": "Wheat"}

# 1. --- DATA EXTRACTION ---
water_results = {b: {cat: 0.0 for cat in FinalCategories} for b in Studyareas}
crit_water_results = {b: {cat: 0.0 for cat in FinalCategories} for b in Studyareas}
n_results = {b: {cat: 0.0 for cat in FinalCategories} for b in Studyareas}
crit_n_results = {b: {cat: 0.0 for cat in FinalCategories} for b in Studyareas}
p_results = {b: {cat: 0.0 for cat in FinalCategories} for b in Studyareas}
crit_p_results = {b: {cat: 0.0 for cat in FinalCategories} for b in Studyareas}

for basin in Studyareas:
    low_runoff_path = os.path.join(data_dir, basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]

    for crop in InputCrops:
        file_path = os.path.join(input_dir, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(file_path): continue
        target_cat = crop_mapper[crop]
        
        with xr.open_dataset(file_path) as ds:
            mask = ds['Basin_mask']; t_m = ds['Total_HA'] > 2500
            mask_not_low_runoff = xr.where(low_runoff.isnull() & (mask == 1), 1, np.nan)

            water_results[basin][target_cat] += (1e-9 *ds['Total_irrigation_amount'].where(t_m) * mask_not_low_runoff * mask).sum().item()
            crit_water_results[basin][target_cat] += (1e-9 * ds['Sus_irrigation_amount'].where(t_m) * mask_not_low_runoff * mask).sum().item()
            n_results[basin][target_cat] += (1e-6 * ds['N_Runoff'].where(t_m) * mask_not_low_runoff * mask).sum().item()
            crit_n_results[basin][target_cat] += (1e-6 * ds['Crit_N_Runoff'].where(t_m) * mask_not_low_runoff * mask).sum().item()
            p_results[basin][target_cat] += (1e-6 * ds['P_Runoff'].where(t_m) * mask_not_low_runoff * mask).sum().item()
            crit_p_results[basin][target_cat] += (1e-6 * ds['Crit_P_Runoff'].where(t_m) * mask_not_low_runoff * mask).sum().item()
            
# 2. --- PLOTTING FUNCTION ---
def plot_1x4_nutrient(data_dict, crit_dict, nutrient_label, filename):
    fig, axes = plt.subplots(1, 4, figsize=(20, 8), sharey=False)
    plt.subplots_adjust(wspace=0.35, top=0.88, bottom=0.2)
    bar_width = 0.4
    
    for i, basin in enumerate(Studyareas):
        ax = axes[i]
        b_actual, b_crit = 0.0, 0.0
        
        for cat in FinalCategories:
            v_act = data_dict[basin][cat]; v_cri = crit_dict[basin][cat]
            ax.bar(0, v_act, bottom=b_actual, width=bar_width, color=crop_colors[cat], edgecolor='white', linewidth=0.5)
            ax.bar(1, v_cri, bottom=b_crit, width=bar_width, color=crop_colors[cat], edgecolor='white', linewidth=0.5)
            b_actual += v_act; b_crit += v_cri
            
        # --- "Nice Number" Axis Logic ---
        max_h = max(b_actual, b_crit)
        if max_h > 0:
            # Determine appropriate step based on magnitude
            if max_h <= 0.5: base = 0.1
            elif max_h <= 1: base = 0.2
            elif max_h <= 2: base = 0.5
            elif max_h <= 5: base = 1
            elif max_h <= 10: base = 2
            elif max_h <= 20: base = 5
            elif max_h <= 50: base = 10
            elif max_h <= 100: base = 20
            elif max_h <= 200: base = 50
            elif max_h <= 500: base = 100
            elif max_h <= 1000: base = 200
            else: base = 500
            
            # Set ymax to the next clean multiple of 'base'
            ymax = np.ceil(max_h / base) * base
            ax.set_ylim(0, ymax)
            
            # Divide ticks nicely based on the same base
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base))
            ax.tick_params(axis='both', labelsize=35)

        # Formatting
        ax.set_xlim(-0.6, 1.6)
        ax.set_xticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
    plt.close()

# 3. --- GENERATE ---
plot_1x4_nutrient(n_results, crit_n_results, "Nitrogen (N)", "Fig4a_N_Boundary_Exceedance_Bar.png")
plot_1x4_nutrient(p_results, crit_p_results, "Phosphorus (P)", "Fig4b_P_Boundary_Exceedance_Bar.png")
plot_1x4_nutrient(water_results, crit_water_results, "Water", "Fig4c_Water_Boundary_Exceedance_Bar.png")

print("Figures saved with clean Y-axis limits.")