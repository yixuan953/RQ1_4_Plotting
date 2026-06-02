# This code is used to plot hte barplots for boudanries for total N and P delivery, agri N and P load, and cropland N and P load in [ktons]

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import os


csv_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Statistics/1_Boundary"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig2"
basins = ["Indus", "LaPlata", "Yangtze", "Rhine"]


crop_colors = {
    'Wheat': "#EBBB7C",
    'Maize': "#B36767",
    'Rice': "#C89CDC",
    'Soybean': "#42949F"
}


def create_boundary_plot(basin, element):
    file_path = os.path.join(csv_dir, f"{basin}_crop-specific_boundaries.csv")
    if not os.path.exists(file_path): return
    
    df = pd.read_csv(file_path, index_col=0)
    col_name = f"{element} [ktons]"
    
    # 1. Values
    total_val = df.loc['All sources', col_name]
    agri_val = df.loc['Agriculture', col_name]
    cropland_val = df.loc['Cropland', col_name]
    
    # 2. Crop Values
    crop_rows = df.reindex(['winterwheat', 'maize', 'mainrice', 'secondrice', 'soybean'])
    crop_rows = crop_rows[col_name].fillna(0)

    crops = {
        'Soybean': crop_rows['soybean'],
        'Rice': crop_rows['mainrice'] + crop_rows['secondrice'],
        'Maize': crop_rows['maize'],
        'Wheat': crop_rows['winterwheat']
    }

    # 3. Plotting
    fig, ax = plt.subplots(figsize=(2.5, 8)) 
    
    pos = 0
    full_width = 0.6 
    edge_style_all_sources = {'edgecolor': '#D7EBF3', 'linewidth': 0.8}
    edge_style_major_crops = {'edgecolor': "#D2D2D2", 'linewidth': 1.2} 
    
    # Draw the main "All sources" bar as the background container
    ax.bar(pos, total_val, color="#D7EBF3", width=full_width, zorder=1, **edge_style_all_sources)
    
    # Draw Agriculture and Cropland as horizontal lines across the bar
    # We calculate the x-min and x-max based on the bar's position and width
    x_min = pos - (full_width / 2)
    x_max = pos + (full_width / 2)
    
    # Define colors for the reference lines
    # Example: Darker blue for Agriculture, Slate/Grey for Cropland
    line_colors = ["#26649E", "#133C63"] 
    
    # Draw Agriculture and Cropland as horizontal lines across the bar
    x_min = pos - (full_width / 2)
    x_max = pos + (full_width / 2)
    
    # ax.hlines accepts a list for 'colors' corresponding to the list of y-values
    ax.hlines([agri_val, cropland_val], x_min, x_max, 
              colors=line_colors, linestyles='solid', linewidth=2.5, zorder=3)
    
    # Major Crops (Stacked at the very bottom)
    bottom_accum = 0
    for name, val in crops.items():
        if val > 0:
            ax.bar(pos, val, bottom=bottom_accum, color=crop_colors[name], 
                   width=full_width, zorder=4, **edge_style_major_crops)
            bottom_accum += val

    # --- Styling ---
    if total_val > 3000: step = 1000
    elif total_val > 1200: step = 500
    elif total_val > 600: step = 200
    elif total_val > 300: step = 100
    elif total_val > 150: step = 50
    elif total_val > 60: step = 20
    elif total_val > 20: step = 10
    else: step = 5    
    
    upper_bound = int(np.ceil(total_val / step)) * step
    ax.set_ylim(0, upper_bound)
    
    # Y-axis ticks
    ticks = np.arange(0, upper_bound + 1, step)
    ax.set_yticks(ticks)
    ax.set_yticklabels([str(int(t)) for t in ticks])
    ax.tick_params(axis='y', labelsize=25, pad=10)
    
    # Spine and Axis Cleanup
    ax.set_xticks([])
    for spine in ['bottom', 'top', 'right']:
        ax.spines[spine].set_visible(False)
    
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['left'].set_position(('outward', 8))
    ax.set_position([0.47, 0.1, 0.5, 0.8])
    
    # Save
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    save_path = os.path.join(output_dir, f"{basin}_{element}_pillar.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"{basin} -- Done!")

for basin in basins:
    create_boundary_plot(basin, "N")
    create_boundary_plot(basin, "P")