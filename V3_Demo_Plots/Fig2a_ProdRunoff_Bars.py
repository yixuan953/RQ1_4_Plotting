# This scipt is used to plot the
# 1 - Total crop production of rainfed & irrigated field [ktons]
# 2 - Total cropland N runoff to surface runoff of rainfed & irrigated field [ktons]
# 3 - Total cropland P runoff to surface runoff of rainfed & irrigated field [ktons]

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# --- Configuration ---
results_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/2_Simulated_Yield_Runoff/1_Current"
fig_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig2_Sim_Yield_Loss/2a_BarPlots"
basins = ["LaPlata", "Rhine",  "Indus", "Yangtze"]
all_crops = ['Wheat', 'Maize', 'Rice', 'Soybean']

crop_colors = {
    'Wheat': "#FFD87D",
    'Maize': '#C96868',
    'Rice': "#BB8ED0",
    'Soybean': '#3A8B95'
}

metrics = [
    ('Production', 'Total_Production_Rainfed', 'Total_Production_Irrigated', '[ktons]'),
    ('N_Runoff', 'Total_N_runoff_Rainfed', 'Total_N_runoff_Irrigated', '[ktons]'),
    ('P_Runoff', 'Total_P_runoff_Rainfed', 'Total_P_runoff_Irrigated', '[ktons]')
]

plt.rcParams.update({'font.size': 18, 'font.family': 'sans-serif'})

def get_smart_xlim(max_val):
    thresholds = [0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000]
    for t in thresholds:
        if max_val <= t: return t
    return np.ceil(max_val / 10000) * 10000

# --- Main Execution ---
for title, rain_row, irri_row, unit in metrics:
    # 1. First pass: Determine how many crops each basin actually has
    basin_data = []
    height_ratios = []
    conv = 1e-9 if title == 'Production' else 1e-6

    for basin in basins:
        file_path = os.path.join(results_dir, f"{basin}_summary_stats.csv")
        if not os.path.exists(file_path):
            basin_data.append(None)
            height_ratios.append(1) # Placeholder height
            continue
        
        df = pd.read_csv(file_path, index_col=0)
        valid_crops_in_basin = []
        
        for crop in all_crops:
            # Check if crop exists in data and has non-zero values
            r_val = df.loc[rain_row, crop] * conv if crop in df.columns and not pd.isna(df.loc[rain_row, crop]) else 0
            i_val = df.loc[irri_row, crop] * conv if crop in df.columns and not pd.isna(df.loc[irri_row, crop]) else 0
            
            if r_val > 0 or i_val > 0:
                valid_crops_in_basin.append({'name': crop, 'rain': r_val, 'irri': i_val})
        
        basin_data.append(valid_crops_in_basin)
        # The height ratio is simply the number of crops (min 1 to avoid errors)
        height_ratios.append(max(len(valid_crops_in_basin), 1))

    # 2. Create figure with dynamic height ratios
    fig, axes = plt.subplots(4, 1, figsize=(8, 14), 
                             gridspec_kw={'height_ratios': height_ratios})
    
    for idx, (ax, data, basin) in enumerate(zip(axes, basin_data, basins)):
        if data is None or len(data) == 0:
            ax.text(0.5, 0.5, f"No crop data for {basin}", ha='center')
            ax.axis('off')
            continue

        y_labels = [d['name'] for d in data]
        rain_vals = [d['rain'] for d in data]
        irri_vals = [d['irri'] for d in data]
        y_pos = np.arange(len(y_labels))

        # Plot Bars
        ax.barh(y_pos, irri_vals, color=[crop_colors[c] for c in y_labels], 
                edgecolor='white', height=0.7, label='Irrigated')
        ax.barh(y_pos, rain_vals, left=irri_vals, color=[crop_colors[c] for c in y_labels], 
                edgecolor='white', height=0.7, alpha=0.4, label='Rainfed')
        
        # --- Styling ---        
        # X-axis logic
        total_vals = [r + i for r, i in zip(rain_vals, irri_vals)]
        ax.set_xlim(0, get_smart_xlim(max(total_vals)))
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.tick_params(axis='x', labelsize=35)
        

        # Y-axis logic
        ax.set_yticks(y_pos)
        ax.set_yticklabels(y_labels, fontsize=35)
        ax.invert_yaxis()

        # Clean frame
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(True)
        ax.spines['left'].set_visible(True)

    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    
    if not os.path.exists(fig_dir): os.makedirs(fig_dir)
    save_path = os.path.join(fig_dir, f"Fig2_{title}_DynamicRows.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

print("Figures generated: Subplot heights adjusted dynamically to crop counts.")
