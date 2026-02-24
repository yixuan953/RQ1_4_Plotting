# This code is used to compare the crop production of 4 scenarios

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

# --- Configuration ---
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
InputCrops = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
FinalCategories = ["Wheat", "Maize", "Rice", "Soybean"]

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
sc1_input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
sc2_input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/2_Sus_irrigation"
sc3_input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
sc4_input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/4_Inc_fert"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs"

cat_colors = {
    "Baseline": "#D97041", 
    "Sustainable Irrigation": "#4A90A2", 
    "Sustainable Irrigation + Reduced Fertilizers": "#7FB069",
    "Sustainable Irrigation + Reduced Fertilizers + Allowable increase": "#C9E743EF"
}
scenarios = list(cat_colors.keys())
crop_mapper = {"winterwheat": "Wheat", "maize": "Maize", "mainrice": "Rice", "secondrice": "Rice", "soybean": "Soybean"}

def lighten_color(color, amount=0.3): # The color will be lighter when amount is higher
    """Lightens the given hex color by mixing it with white."""
    try:
        c = mcolors.cnames[color]
    except:
        c = color
    c = mcolors.to_rgb(c)
    return mcolors.to_hex([1 - amount * (1 - x) for x in c])

# --- Data Processing and Plotting ---
for basin in Studyareas:
    # Dictionary to store results: {Crop: {Scenario: {'rainfed': val, 'irrigated': val}}}
    plot_data = {cat: {sc: {'rainfed': 0.0, 'irrigated': 0.0} for sc in scenarios} for cat in FinalCategories}

    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)

    for crop in InputCrops:
        sc1_nc = os.path.join(sc1_input_dir, f"{basin}_{crop}_summary_baseline.nc")
        sc2_nc = os.path.join(sc2_input_dir, f"{basin}_{crop}_summary.nc")
        sc3_nc = os.path.join(sc3_input_dir, f"{basin}_{crop}_summary.nc")
        sc4_nc = os.path.join(sc4_input_dir, f"{basin}_{crop}_summary.nc")

        if not all(os.path.exists(f) for f in [sc1_nc, sc2_nc, sc3_nc, sc4_nc]):
            print(f"Skipping {basin} - {crop}: Files missing.")
            continue

        target_cat = crop_mapper[crop]
        
        for sc_idx, (sc_nc, sc_label) in enumerate(zip([sc1_nc, sc2_nc, sc3_nc, sc4_nc], scenarios)):
            ds = xr.open_dataset(sc_nc)
            mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)

            # Rainfed calculation
            prod_rf = 0.000001 * (ds["Avg_Yield_Rainfed"] * ds["Rainfed_HA"] * mask * mask_not_low_runoff).sum(skipna=True).values
            # Irrigated calculation (Fixed HA variable from your snippet)
            prod_ir = 0.000001 * (ds["Avg_Yield_Irrigated"] * ds["Irrigated_HA"] * mask * mask_not_low_runoff).sum(skipna=True).values
            
            plot_data[target_cat][sc_label]['rainfed'] += float(prod_rf)
            plot_data[target_cat][sc_label]['irrigated'] += float(prod_ir)
            ds.close()
            
    # --- Print Results to Console ---
    print(f"\n{'='*60}")
    print(f" BASIN: {basin}")
    print(f"{'Crop':<10} | {'Scenario':<30} | {'Rainfed':>10} | {'Irrig.':>10}")
    print(f"{'-'*60}")
    for cat in FinalCategories:
        for sc in scenarios:
            rf = plot_data[cat][sc]['rainfed']
            ir = plot_data[cat][sc]['irrigated']
            print(f"{cat:<10} | {sc[:30]:<30} | {rf:10.8f} | {ir:10.8f}")
    print(f"{'='*60}")


    # --- Plotting Section ---
    fig, ax = plt.subplots(figsize=(14, 6))
    
    n_cats = len(FinalCategories)
    n_scens = len(scenarios)
    bar_width = 0.2
    group_gap = 0.1  # Gap between different crop types
    
    # Calculate x positions
    indices = np.arange(n_cats)
    
    for i, sc_label in enumerate(scenarios):
        # Position for each scenario bar within the group
        pos = indices + (i - (n_scens - 1)/2) * bar_width
        
        rf_vals = [plot_data[cat][sc_label]['rainfed'] for cat in FinalCategories]
        ir_vals = [plot_data[cat][sc_label]['irrigated'] for cat in FinalCategories]
        
        # Color definitions
        base_color = cat_colors[sc_label]
        light_color = lighten_color(base_color, 0.4) # Slightly lighter for rainfed
        
        # Plot stacked bars
        # Bottom: Irrigated (stacked on rainfed)
        ax.bar(pos, ir_vals, bar_width, label=f"{sc_label} (Irrigated)" if i==0 else "", 
               color=base_color, edgecolor='white', linewidth=0.5)
        # Top: Rainfed 
        ax.bar(pos, rf_vals, bar_width, bottom=ir_vals, label=f"{sc_label} (Rainfed)" if i==0 else "", 
               color=light_color, edgecolor='white', linewidth=0.5)

    # Styling
    ax.set_title(f"Crop Production - {basin} Basin", fontsize=14, fontweight='bold')
    ax.set_ylabel("ktons", fontsize=12)
    ax.set_xticks(indices)
    ax.set_xticklabels(FinalCategories, fontsize=12)
    
    # Create custom legend to avoid duplicate labels
    legend_elements = []
    for sc_label in scenarios:
        legend_elements.append(Line2D([0], [0], color=cat_colors[sc_label], lw=6, label=sc_label))
    
    # Add a note about the shading in legend or title
    ax.legend(handles=legend_elements, title="Scenarios", loc='upper left', bbox_to_anchor=(1, 1))
    ax.text(1.02, 0.5, "Lighter shade: Rainfed\nDarker shade: Irrigated", 
            transform=ax.transAxes, verticalalignment='center')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basin}_Production_Comparison.png"), dpi=300)
    plt.close()

print("Processing complete. Plots saved to output directory.")