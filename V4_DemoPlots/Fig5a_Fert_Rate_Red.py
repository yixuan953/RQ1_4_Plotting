# This script is used to plot the changes of fertilization rate [kg/ha]

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import os

# Paths
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Statistics/4_Fert_Red"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig5"
Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]

def plot_butterfly_fertilization(study_areas, input_path, output_path):
    # 1. Load data
    all_dfs = []
    for basin in study_areas:
        csv_file = os.path.join(input_path, f"{basin}_Fert_Final_Summary.csv")
        if os.path.exists(csv_file):
            temp_df = pd.read_csv(csv_file)
            temp_df['Basin'] = basin
            all_dfs.append(temp_df)
    
    if not all_dfs: return
    
    # 2. Dynamic Height
    total_crops = sum(len(d) for d in all_dfs)
    inch_per_crop = 1.8
    fig_height = (total_crops * inch_per_crop) + 4
    
    # Create two subplots side-by-side with shared Y axis
    # width_ratios can be adjusted if N and P need different relative widths
    fig, (ax_n, ax_p) = plt.subplots(1, 2, figsize=(14, fig_height), sharey=True)
    plt.subplots_adjust(wspace=0.4) # Adjust space for crop labels in the middle
    plt.rcParams['font.family'] = 'DEJAVU SANS'
    plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans']

    # Colors and Settings
    bar_height = 0.5
    offset = bar_height / 2
    c_inorg, c_manure, p_inorg, p_manure =  "#C1864B", "#C1864B",  "#A3CD7C", "#A3CD7C"
    
    current_y_pos = 0
    ytick_positions = []
    ytick_labels = []

    # 3. Plotting Loop
    for df in all_dfs:
        basin_name = df['Basin'].iloc[0]
        num_crops = len(df)
        y_indices = np.arange(current_y_pos, current_y_pos + num_crops)
        
        for nutrient, ax in [('N', ax_n), ('P', ax_p)]:
            # Column names
            c_in = f'Current_Rate_{nutrient}_Inorg'
            c_ma = f'Current_Rate_{nutrient}_Manure'
            p_in = f'PostRed_Rate_{nutrient}_Inorg'
            p_ma = f'PostRed_Rate_{nutrient}_Manure'
            
            # --- Plot Current (Top bar) ---
            # Inorganic
            ax.barh(y_indices - offset, df[c_in], bar_height, 
                    color=c_inorg, edgecolor='grey')
            # Manure (Added hatch and same color)
            ax.barh(y_indices - offset, df[c_ma], bar_height, left=df[c_in], 
                    color=c_manure, edgecolor='white', hatch='///') 
            
            # --- Plot Post-Red (Bottom bar) ---
            # Inorganic
            ax.barh(y_indices + offset, df[p_in], bar_height, 
                    color=p_inorg, edgecolor='grey')
            # Manure (Added hatch and same color)
            ax.barh(y_indices + offset, df[p_ma], bar_height, left=df[p_in], 
                    color=p_manure, edgecolor='white', hatch='///')
            
            # Percentage Annotations
            for i in range(num_crops):
                total_curr = df.iloc[i][c_in] + df.iloc[i][c_ma]
                total_post = df.iloc[i][p_in] + df.iloc[i][p_ma]
                if total_curr > 0:
                    reduction = ((total_curr - total_post) / total_curr) * 100
                    # For N (left side), text moves further left; for P (right side), text moves further right
                    text_x = total_post + 3
                    ax.text(text_x, y_indices[i] + offset, f'-{reduction:.0f}%', 
                            va='center', fontsize=30, color="#2A2929FF",
                            ha='right' if nutrient == 'N' else 'left')

        # Update tick tracking
        ytick_positions.extend(y_indices)
        ytick_labels.extend(df['Crop_Group'])

        current_y_pos += num_crops + 2

    # --- Formatting N (Left Side) ---
    ax_n.set_xlim(0, 300)
    ax_n.invert_xaxis() # This creates the middle-to-left increase
    ax_n.xaxis.tick_top()
    ax_n.xaxis.set_label_position('top')
    ax_n.set_xlabel('N Fertilization Rate\n[$kg ha^{-1} yr^{-1}$]', fontsize=30, labelpad=20)
    ax_n.tick_params(axis='x', labelsize=30)
    ax_n.spines['left'].set_visible(False)
    ax_n.spines['bottom'].set_visible(False)
    ax_n.spines['right'].set_visible(True) # Border in the middle

    # --- Formatting P (Right Side) ---
    ax_p.set_xlim(0, 60)
    ax_p.xaxis.tick_top()
    ax_p.xaxis.set_label_position('top')
    ax_p.set_xlabel('P Fertilization Rate\n[$kg ha^{-1} yr^{-1}$]', fontsize=30, labelpad=20)
    ax_p.tick_params(axis='x', labelsize=30)
    ax_p.spines['right'].set_visible(False)
    ax_p.spines['bottom'].set_visible(False)

    # --- Middle Y-Axis (Crop Labels) ---
    # We remove labels from subplots and place them in the middle
    ax_n.set_yticks(ytick_positions)
    ax_n.set_yticklabels([]) # Hide default
    ax_p.set_yticklabels([]) # Hide default shared
    ax.tick_params(axis='both', labelsize=30)
    
    for y_pos, label in zip(ytick_positions, ytick_labels):
        # Place text in the white space between subplots
        ax_p.text(-0.2, y_pos, label, transform=ax_p.get_yaxis_transform(),
                  ha='center', va='center', fontsize=30)

    ax_n.invert_yaxis()
    
    save_path = os.path.join(output_path, "Butterfly_N_P_FertRate.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print("Butterfly plot generated.")

# Run
plot_butterfly_fertilization(Studyareas, input_dir, output_dir)
