# This scipt is used to plot the
# 1 - Total crop production of rainfed & irrigated field [ktons]
# 2 - Total cropland N runoff to surface runoff of rainfed & irrigated field [ktons]
# 3 - Total cropland P runoff to surface runoff of rainfed & irrigated field [ktons]

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator


# --- Configuration ---
results_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Statistics/2_Simulated_Yield_Runoff"
fig_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig3"
basins = ["LaPlata", "Rhine",  "Indus", "Yangtze"]
all_crops = ['Wheat', 'Maize', 'Rice', 'Soybean']

crop_colors = {
    'Wheat': "#E8AE62",
    'Maize': "#B55454",
    'Rice': "#C48BDE",
    'Soybean': "#3B939F"
}

metrics = [
    ('Production', 'Total_Production_Rainfed', 'Total_Production_Irrigated', '[ktons]'),
    ('N_Runoff', 'Total_N_runoff_Rainfed', 'Total_N_runoff_Irrigated', '[ktons]'),
    ('P_Runoff', 'Total_P_runoff_Rainfed', 'Total_P_runoff_Irrigated', '[ktons]')
]

plt.rcParams.update({'font.family': 'DejaVu Sans'})


def get_smart_xlim(max_val):
    thresholds = [0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20, 25, 40, 60, 80, 100, 200, 400, 600, 800, 1000]
    for t in thresholds:
        if max_val <= t: return t
    return np.ceil(max_val / 1000) * 1000

for title, rain_row, irri_row, unit in metrics:
    # 1. Data Prep and Width Calculation
    basin_data = []
    num_crops_per_basin = []
    conv = 1e-9 if title == 'Production' else 1e-6

    for basin in basins:
        file_path = os.path.join(results_dir, f"{basin}_summary_stats.csv")
        if not os.path.exists(file_path):
            basin_data.append(None)
            num_crops_per_basin.append(1) # Placeholder for empty plot
            continue
        
        df = pd.read_csv(file_path, index_col=0)
        valid_crops_in_basin = []
        for crop in all_crops:
            r_val = df.loc[rain_row, crop] * conv if crop in df.columns and not pd.isna(df.loc[rain_row, crop]) else 0
            i_val = df.loc[irri_row, crop] * conv if crop in df.columns and not pd.isna(df.loc[irri_row, crop]) else 0
            if r_val > 0 or i_val > 0:
                valid_crops_in_basin.append({'name': crop, 'rain': r_val, 'irri': i_val})
        
        basin_data.append(valid_crops_in_basin)
        # Ratio is based on crop count; this makes 2 bars look 1/5th the size of 10 bars
        num_crops_per_basin.append(max(len(valid_crops_in_basin), 1))

    # 2. Create figure with Width Ratios
    # Adjust figsize width (22) if you have many crops total
    fig = plt.figure(figsize=(22, 11))
    gs = gridspec.GridSpec(1, 4, width_ratios=num_crops_per_basin)
    axes = [fig.add_subplot(gs[i]) for i in range(4)]

    for idx, (ax, data, basin) in enumerate(zip(axes, basin_data, basins)):
        if data is None or len(data) == 0:
            ax.axis('off')
            continue

        y_labels = [d['name'] for d in data]
        rain_vals = np.array([d['rain'] for d in data])
        irri_vals = np.array([d['irri'] for d in data])
        x_pos = np.arange(len(y_labels))

        # Vertical bars touching (width=1.0)
        ax.bar(x_pos, irri_vals, color=[crop_colors[c] for c in y_labels], 
               edgecolor='grey', width=1.0, linewidth=0.8, zorder=3)
        ax.bar(x_pos, rain_vals, bottom=irri_vals, color=[crop_colors[c] for c in y_labels], 
               edgecolor='grey', width=1.0, linewidth=0.8, alpha=0.3, zorder=3)
        
        # --- Y-Axis Styling ---
        max_total = max(rain_vals + irri_vals)
        upper_limit = get_smart_xlim(max_total)
        ax.set_ylim(0, upper_limit)
        
        # Limit ticks to max 5
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        ax.tick_params(axis='y', labelsize=50)

        # --- X-Axis Styling ---
        # Locking limits tightly to the bars is what fixes the width consistency
        ax.set_xlim(-0.5, len(y_labels) - 0.5)
        ax.set_xticks([])

        # --- Frame ---
        for spine in ['top', 'right']: ax.spines[spine].set_visible(True)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

    # Prevent y-labels from causing horizontal shifts between subplots
    fig.align_ylabels(axes)
    plt.tight_layout() 
    plt.subplots_adjust(wspace=0.8)

    save_path = os.path.join(fig_dir, f"Fig3_{title}_Fixed_Barplots.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

print("Figures generated: Subplot heights adjusted dynamically to crop counts.")
