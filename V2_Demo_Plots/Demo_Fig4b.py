# This code is used to plot the changes in N, P runoff for different sceanrios (Bar plot)
from matplotlib.lines import Line2D
import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


# --- Configuration ---
# Rice
# Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
# InputCrops = ["mainrice", "secondrice"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
# output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4b/Rice"

# Soybean
Studyareas = ["LaPlata", "Yangtze"]
InputCrops = ["soybean"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4b/Soybean"

# # # Winterwheat
# Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
# InputCrops = ["winterwheat"]
# output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4b/Wheat"

# # Maize
# Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
# InputCrops = ["maize"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
# output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4b/Maize"

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
sc_dirs = [
    "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline",
    "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/2_Sus_irrigation",
    "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert",
    "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/4_Inc_fert"
]

os.makedirs(output_dir, exist_ok=True)

cat_colors = ["#EB8B39", "#65A2D2", "#72B285", "#D4D98A"]
sc_labels = ["Baseline", "Sus Irri", "Sus Irri + Red Fert", "Sus Irri + Red Fert + Allowable Inc"]
    
# --- Updated Plotting Function: Slim, No Title, No X-Ticks, with Legend ---
def save_separate_waterfall(basin, values, thresh, var_name, color_list, sc_labels, output_dir):
    # Set a slim, portrait-oriented figsize
    fig, ax = plt.subplots(figsize=(4, 8))
    
    x = np.arange(len(values))
    
    # 1. Plotting the Bars (Baseline and Step Changes)
    for i in range(len(values)):
        if i == 0:
            # Step bars: removed edgecolor by setting linewidth=0
            ax.bar(x[i], values[i], color=color_list[i], lw=0, width=0.5, zorder=3)
        else:
            prev_val = values[i-1]
            curr_val = values[i]
            bottom = min(prev_val, curr_val)
            height = abs(prev_val - curr_val)
            
            ax.bar(x[i], height, bottom=bottom, color=color_list[i], 
                            lw=0, width=0.5, zorder=3)
            
            # Connecting line
            ax.plot([x[i-1], x[i]], [prev_val, prev_val], color='gray', 
                    linestyle=':', linewidth=1.5, zorder=2)

    # 2. Critical Threshold Line
    ax.axhline(y=thresh, color="#f06a5b", linestyle='--', linewidth=2.5, zorder=4)

    # 3. Styling: Remove top, right, and bottom borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Optional: remove bottom spine if you want it truly floating
    ax.spines['bottom'].set_visible(False) 

    ax.tick_params(axis='y', labelsize=15)
    ax.set_ylabel(f'Annual Production (ktons)', fontsize=15)
    
    # --- REMOVE X-TICKS AND LABELS ---
    ax.set_xticks([]) 
    ax.set_xticklabels([])

    # 4. Create Custom Legend
    # We create "proxies" (rectangles for bars and a line for the threshold)
    legend_elements = [
        plt.Rectangle((0,0),1,1, color=color_list[i], ec='k', label=sc_labels[i]) 
        for i in range(len(sc_labels))
    ]
    legend_elements.append(Line2D([0], [0], color='#e74c3c', lw=2.5, ls='--', label='Critical Value'))

    # Place legend below the plot
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
              ncol=1, frameon=False, fontsize=15)

    # Headroom for the y-axis
    ax.set_ylim(0, max(values.max(), thresh) * 1.15)
    
    # Light grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, zorder=0)

    # Adjust layout to make room for the legend at the bottom
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    
    # Save files separately
    file_name = f"{basin}_Waterfall_{var_name}.png"
    plt.savefig(os.path.join(output_dir, file_name), dpi=300, transparent=True)
    plt.close()


for basin in Studyareas:
    print(f"Processing Basin: {basin}")
    
    # Load low runoff mask
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        mask_not_low_runoff = xr.where(ds_lr["Low_Runoff"].isnull(), 1, np.nan)

    # --- Define threshold percentage ---
    threshold_lookup = {
        "Yangtze": 0.75,
        "Indus": 0.75,
        "Rhine": 0.65,
        "LaPlata": 0.20
    }
    threshold_pct = threshold_lookup.get(basin, 1.0) # Default to 100% if not found

    # Initialize storage for basin-wide sums across all crops for each scenario
    total_prod_by_scenario = np.zeros(len(sc_dirs))

    for sc_idx, sc_input_dir in enumerate(sc_dirs):
        for crop in InputCrops:
            nc_path = os.path.join(sc_input_dir, f"{basin}_{crop}_summary.nc")
            if not os.path.exists(nc_path): 
                continue
            
            with xr.open_dataset(nc_path) as ds:
                mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
                combined_mask = mask * mask_not_low_runoff
                
                # Calculate Production (Sum of Irrigated and Rainfed)
                prod_map = (ds["Irrigated_HA"] * ds["Avg_Yield_Irrigated"]) + (ds["Rainfed_HA"] * ds["Avg_Yield_Rainfed"])
                
                # Convert to Million Tons (or your preferred unit) and sum
                crop_prod_sum = (prod_map * combined_mask).sum().item() * 1e-6
                total_prod_by_scenario[sc_idx] += crop_prod_sum

    # --- Calculate the Final Basin Target ---
    # The target is based on the total production of ALL crops in Scenario 0 (Baseline)
    basin_target = total_prod_by_scenario[0] * threshold_pct

    print(f"Basin: {basin} | Baseline: {total_prod_by_scenario[0]:.2f} | Target: {basin_target:.2f}")

    # --- Alert Logic ---
    for sc_idx, total in enumerate(total_prod_by_scenario):
        if sc_idx > 0 and total < basin_target:
            print(f"⚠️ Alert: Scenario {sc_idx} total ({total:.2f}) is below target!")

    # Save the waterfall plot using the aggregate target
    save_separate_waterfall(basin, total_prod_by_scenario, basin_target, 'CropProd', cat_colors, sc_labels, output_dir)