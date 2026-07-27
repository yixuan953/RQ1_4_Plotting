import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches


# --- Configuration (Maintained) ---
Studyareas = ["LaPlata", "Rhine", "Indus"] #"Yangtze"]
InputCrops = ["soybean", "mainrice", "secondrice", "maize", "winterwheat"]
FinalCategories = ["Soybean", "Rice", "Maize", "Wheat"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig2"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"

# REVISED: Integrated pattern mapping & background colors as requested
crop_hatches = {
    'Wheat': '||',     # Less dense vertical lines
    'Maize': 'xx',    # Less dense diagonal grid
    'Rice': '..',      # Less dense dots
    'Soybean': '\\\\'  # Less dense backward slashes
}
crop_mapper = {"soybean": "Soybean",  "mainrice": "Rice", "secondrice": "Rice", "maize": "Maize", "winterwheat": "Wheat"}

COLOR_CURRENT = "#D3D3D3"  # Light Grey for current bars
COLOR_CRIT = "#FFFFFF"     # Light Blue for critical bars

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
    fig, axes = plt.subplots(1, 4, figsize=(7.5, 6), sharey=False)
    plt.subplots_adjust(wspace=0.41, top=0.88, bottom=0.2)
    bar_width = 0.3
    
    for i, basin in enumerate(Studyareas):
        ax = axes[i]
        b_actual, b_crit = 0.0, 0.0
        
        for cat in FinalCategories:
            v_act = data_dict[basin][cat]; v_cri = crit_dict[basin][cat]
            hatch_pattern = crop_hatches[cat]
            
            # REVISED LOGIC: Uniform base color fills combined with clean black hatch patterns for stacking visibility
            ax.bar(0 - bar_width/2, v_act, bottom=b_actual, width=bar_width, 
                   color=COLOR_CURRENT, edgecolor="#3B3B3B", linewidth=2.0, hatch=hatch_pattern)
            
            ax.bar(0 + bar_width/2, v_cri, bottom=b_crit, width=bar_width, 
                   color=COLOR_CRIT, edgecolor='#3B3B3B', linewidth=2.0, hatch=hatch_pattern)
            
            b_actual += v_act; b_crit += v_cri
            
        # --- "Nice Number" Axis Logic (Allows Clean Decimals for Totals < 2) ---
        max_h = max(b_actual, b_crit)
        if max_h > 0:
            if max_h <= 0.03: base = 0.01      # Ticks: 0, 0.01, 0.02, 0.03
            elif max_h <= 0.06: base = 0.02    # Ticks: 0, 0.02, 0.04, 0.06
            elif max_h <= 0.15: base = 0.05    # Ticks: 0, 0.05, 0.10, 0.15
            elif max_h <= 0.3: base = 0.1      # Ticks: 0, 0.1, 0.2, 0.3
            elif max_h <= 0.6: base = 0.2      # Ticks: 0, 0.2, 0.4, 0.6
            elif max_h <= 1.5: base = 0.5      # Ticks: 0, 0.5, 1.0, 1.5
            elif max_h <= 3: base = 1          # Ticks: 0, 1, 2, 3
            elif max_h <= 6: base = 2          # Ticks: 0, 2, 4, 6
            elif max_h <= 15: base = 5         # Ticks: 0, 5, 10, 15
            elif max_h <= 30: base = 10        # Ticks: 0, 10, 20, 30
            elif max_h <= 60: base = 20        # Ticks: 0, 20, 40, 60
            elif max_h <= 150: base = 50       # Ticks: 0, 50, 100, 150
            elif max_h <= 300: base = 100      # Ticks: 0, 100, 200, 300
            elif max_h <= 600: base = 200      # Ticks: 0, 200, 400, 600
            elif max_h <= 1500: base = 500     # Ticks: 0, 500, 1000, 1500
            else: base = 1000                  # Ticks: 0, 1000, 2000, 3000
            
            ymax = base * 3
            # In case data outstrips initial bounds, upscale cleanly by increments of 2x or 5x
            while ymax < max_h:
                base *= 2
                ymax = base * 3
                
            ax.set_ylim(0, ymax)
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base))
            ax.tick_params(axis='both', labelsize=20)
            
            # Formatter '%g' drops extraneous trailing floating decimals (.0) elegantly
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
            
            for label in ax.get_yticklabels():
                label.set_fontname('sans-serif')
        
        # Formatting (Maintained)
        ax.set_xlim(-0.6, 0.6)
        ax.set_xticks([]) 
        ax.set_xticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)


    legend_elements = [
            # Part 1: Crop Hatch Types (using a white background to clearly showcase the patterns)
        mpatches.Patch(facecolor='white', edgecolor='#4A4A4A', linewidth=2.0,hatch=crop_hatches['Soybean'], label='Soybean'),
        mpatches.Patch(facecolor='white', edgecolor='#4A4A4A', linewidth=2.0,hatch=crop_hatches['Rice'], label='Rice'),
        mpatches.Patch(facecolor='white', edgecolor='#4A4A4A', linewidth=2.0,hatch=crop_hatches['Maize'], label='Maize'),
        mpatches.Patch(facecolor='white', edgecolor='#4A4A4A', linewidth=2.0,hatch=crop_hatches['Wheat'], label='Wheat'),
            
        # Blank spacer patch to split the legend sections neatly
        mpatches.Patch(color='none', label=''), 
            
           # Part 2: Evaluation Scenarios
        mpatches.Patch(facecolor=COLOR_CURRENT, edgecolor='#4A4A4A',linewidth=2.0, label='Current'),
        mpatches.Patch(facecolor=COLOR_CRIT, edgecolor='#4A4A4A', linewidth=2.0,label='Critical')]

        # Render the global master legend on the figure level
    fig.legend(
        handles=legend_elements, 
        loc='center right', 
        bbox_to_anchor=(0.99, 0.5), # Tweaked slightly right to accommodate larger text
        frameon=True, 
        edgecolor="#FFFFFF", 
        facecolor='white',
        labelspacing=1.2,       # Increased vertical spacing between rows
        handlelength=5,       # Made the sample pattern boxes wider
        handleheight=5        # Made the sample pattern boxes taller
    )
    
    os.makedirs(output_dir, exist_ok=True) 
    plt.savefig(os.path.join(output_dir, filename), dpi=300) # ,bbox_inches='tight')
    plt.close()

# 3. --- GENERATE ---
plot_1x4_nutrient(n_results, crit_n_results, "Nitrogen (N)", "Fig4a_N_Boundary_Exceedance_Bar.png")
plot_1x4_nutrient(p_results, crit_p_results, "Phosphorus (P)", "Fig4b_P_Boundary_Exceedance_Bar.png")
plot_1x4_nutrient(water_results, crit_water_results, "Water", "Fig4c_Water_Boundary_Exceedance_Bar.png")

print("Figures successfully saved with matching side-by-side patterned segments.")