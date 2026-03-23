# This code is used to plot hte barplots for boudanries for total N and P load, agri N and P load, and cropland N and P load in [ktons]

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import os


csv_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/1_Boundary_load"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig1_Boundary/1a_Bars"
basins = ["Indus", "LaPlata", "Yangtze", "Rhine"]


crop_colors = {
    'Wheat': "#FADFA1",
    'Maize': '#C96868',
    'Rice': "#BB8ED0",
    'Soybean': '#3A8B95'
}

def create_boundary_plot(basin, element):
    file_path = os.path.join(csv_dir, f"{basin}_crit_load_sum_updated.csv")
    if not os.path.exists(file_path): return
    
    df = pd.read_csv(file_path, index_col=0)
    col_name = f"{element} [ktons]"
    
    # 1. Values
    total_val = df.loc['All sources', col_name]
    agri_val = df.loc['Agriculture', col_name]
    cropland_val = df.loc['Cropland', col_name]
    
    # 2. Crop Values
    wheat = df.loc['winterwheat', col_name] if 'winterwheat' in df.index else 0
    maize = df.loc['maize', col_name] if 'maize' in df.index else 0
    rice = (df.loc['mainrice', col_name] if 'mainrice' in df.index else 0) + \
           (df.loc['secondrice', col_name] if 'secondrice' in df.index else 0)
    soybean = df.loc['soybean', col_name] if 'soybean' in df.index else 0
    
    major_crop_sum = wheat + maize + rice + soybean
    crops = {'Wheat': wheat, 'Maize': maize, 'Rice': rice, 'Soybean': soybean}
    
    # 3. Proportions
    prop_agri_total = (agri_val / total_val) * 100
    prop_allcrop = (cropland_val / agri_val) * 100
    prop_majorcrops = (major_crop_sum / cropland_val) * 100

    # 4. Plotting
    fig, ax = plt.subplots(figsize=(10, 3))
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans']
    plt.tick_params(axis='x', labelsize=16)
    y_labels = ['Major crops', 'Cropland', 'Agriculture', 'All sources']
    
    # Bars
    ax.barh(3, total_val, color="#547792", edgecolor= "white", height=0.6, linewidth=2) # All sources
    ax.barh(2, agri_val, color="#AACDDC", edgecolor= "white", height=0.6, linewidth=2) # Agriculture
    ax.barh(1, cropland_val, color="#D6DAC8", edgecolor="white", height=0.6, linewidth=2) # Cropland
    
    left_start = 0
    for name, val in crops.items():
        ax.barh(0, val, left=left_start, color=crop_colors[name], 
                edgecolor="#f8f8f8f9", height=0.6)
        left_start += val
        
    # --- Axes and Frame Styling ---
    
    # Move X-axis to the top
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    
    # Remove Right and Bottom spines (frame), keep Top and Left (axes)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(True)
    ax.spines['left'].set_visible(True)

    # Create a clean xlim
    if total_val > 3000:
        step = 1000
    elif total_val > 1200:
        step = 500
    elif total_val > 600:
        step = 200
    elif total_val > 300:
        step = 100
    elif total_val > 150:
        step = 50
    elif total_val > 60:
        step = 20
    elif total_val > 20:
        step = 10
    elif total_val > 15:
        step = 5
    else: 
        step = 2    
    upper_bound = int(np.ceil(total_val / step)) * step
    
    # 2. Apply the limit
    ax.set_xlim(0, upper_bound)
    
    # 3. Create ticks at every 'step' interval
    # We use +1 to ensure the upper_bound itself is included in the ticks
    ticks = np.arange(0, upper_bound + 1, step)
    ax.set_xticks(ticks)

    # Adjust labels and font sizes
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(y_labels, fontsize=18)
    ax.set_xlabel(f'Boundaries for {element} delivery to surface water [ktons]', 
                  fontsize=18, labelpad=20)
    
    # Text Annotations
    ax.text(agri_val, 2, f'  {prop_agri_total:.0f}% of all sources', 
            va='center', fontsize=18)
    ax.text(cropland_val, 1, f'  {prop_allcrop:.0f}% of agriculture', 
            va='center', fontsize=18)
    ax.text(major_crop_sum, 0, f'  {prop_majorcrops:.0f}% of all crops', 
            va='center', fontsize=18)
        
    # Optional: ensure the axes lines match the tick color
    ax.spines['top'].set_color('#333333')
    ax.spines['left'].set_color('#333333')

    plt.tight_layout()
    
    # Save
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    save_path = os.path.join(output_dir, f"{basin}_{element}_boundary_load.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"{basin} -- Done!")

for basin in basins:
    create_boundary_plot(basin, "N")
    create_boundary_plot(basin, "P")