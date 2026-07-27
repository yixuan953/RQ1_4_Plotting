# This code is used to plot the crop production under different scenarios

import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- 1. Configuration ---
# Paths & Fractions (Adjust these paths if needed)
input_csv = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/2_NPYield_Change/RespectBoundaries_Yield_Runoff_Summary.csv"
current_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/1_Boundary" 
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig4"
os.makedirs(output_dir, exist_ok=True)

# N or P
# element = "N"
# dissolved_frac = 0.7  # Adjust if your actual fraction is less than 1.0
# var_name = "N [ktons]"
# var_name2 = "N_Runoff_ktons"

element = "P"
dissolved_frac = 0.25  # Adjust if your actual fraction is less than 1.0
var_name = "P [ktons]"
var_name2 = "P_Runoff_ktons"

Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
FinalCategories = ["Soybean", "Rice", "Maize", "Wheat"]

# Mapping required to match summary crops to boundary file keys
CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}

# Explicit scenario order (will be plotted from left to right on X axis)
scenario_order = [
    "Current",
    "Water",
    "Water+N",
    "Water+P",  
    "Water+NP",
]

crop_hatches = {
    'Wheat': '|',     
    'Maize': 'x',    
    'Rice': '.',      
    'Soybean': '\\'    
}

scenario_face_colors = {
    "Water+NP": "#02B585",  # green
    "Water+P" : "#E8A921",  # orange
    "Water+N" : "#D970AA",  # purple-pink
    "Water"   : "#56B4E9",  # sky blue
    "Current" : "#999999"   # gray
}

# --- 2. Load and Clean Data ---
df = pd.read_csv(input_csv)

# Filter for only target crops
df = df[df['Crop'].isin(FinalCategories)]

# --- 3. Loop and Plot for Each Basin ---
for basin in Studyareas:
    print(f"Processing basin: {basin}...")
    
    df_basin = df[df['Basin'] == basin].copy()
    if df_basin.empty:
        print(f"No data found for basin: {basin}, skipping...")
        continue
        
    # --- 4. Boundary Loading ---
    available_summary_crops = df_basin['Crop'].unique()
    basin_boundary_crops = []
    for s_crop in available_summary_crops:
        if s_crop in CROP_MAP:
            basin_boundary_crops.extend(CROP_MAP[s_crop])

    current_n_bound = None
    b_file_path = os.path.join(current_boundary_dir, f"{basin}_crop-specific_boundaries.csv")   

    if os.path.exists(b_file_path):
        df_b = pd.read_csv(b_file_path, index_col=0)
        df_b.index = df_b.index.astype(str).str.strip()
        
        active_rows = [row for row in basin_boundary_crops if row in df_b.index]
        if active_rows:
            current_n_bound = dissolved_frac * float(df_b.loc[active_rows, var_name].sum())
    else:
        print(f"  [Warning] Skipping boundary for {basin}: Boundary file missing at {b_file_path}")

    # Pivot the data
    pivot_df = df_basin.pivot(index='Scenario', columns='Crop', values=var_name2)
    
    # Reindex to maintain logical sequence
    pivot_df = pivot_df.reindex(scenario_order)
    pivot_df = pivot_df.reindex(columns=FinalCategories)

    # --- Plotting ---
    # Fix the overall figure canvas size across ALL execution loops
    fig_width = 6.0
    fig_height = 8.0
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # CRITICAL ADJUSTMENT: Force absolute frame box parameters 
    # format: [left_margin, bottom_margin, frame_width, frame_height] as percentage of canvas (0 to 1)
    # This locks the frame to identical physical sizes regardless of label strings.
    ax = fig.add_axes([0.25, 0.1, 0.7, 0.85])
    
    x_indices = np.arange(len(pivot_df.index))
    
    bottom_y = None 
    for crop in FinalCategories:
        if crop in pivot_df.columns:
            lengths = pivot_df[crop]
            
            ax.bar(
                x_indices, 
                lengths, 
                bottom=bottom_y,  
                label=crop, 
                hatch=crop_hatches[crop],
                facecolor=[scenario_face_colors.get(s, 'lightgray') for s in pivot_df.index],
                edgecolor='black',
                linewidth=0.8,
                width=1.0,        # Keeps bars touching flush edge-to-edge
                alpha=0.9,
                zorder=3          
            )
            
            if bottom_y is None:
                bottom_y = lengths.fillna(0)
            else:
                bottom_y += lengths.fillna(0)

    # --- 5. Add Horizontal Boundary Line ---
    if current_n_bound is not None:
        ax.axhline(
            y=current_n_bound, 
            color='red', 
            linestyle='-', 
            linewidth=2.0, 
            label='N Boundary',
            zorder=4  
        )
        print(f"  Added N boundary line at: {current_n_bound} ktons")

    # Force the data frame boundary exactly to the edges of the continuous bars
    ax.set_xlim(-0.5, len(x_indices) - 0.5)
    ax.set_xmargin(0)

    # Clean Up Axes Visuals
    ax.spines['left'].set_color("#000000")
    ax.spines['bottom'].set_color("#000000")
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)
    
    # Hide X labels, keep centered tick placements
    ax.set_xticks(x_indices)
    ax.set_xticklabels([])
    
    ax.tick_params(axis='y', labelsize=35)
    ax.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)  
    ax.set_axisbelow(True)

    # DO NOT use plt.tight_layout() here, as it overrides the absolute positioning set by add_axes
    
    # Save image (using explicit coordinates, padding changes are managed by the margins above)
    out_path = os.path.join(output_dir, f"{basin}_Total_{element}_Runoff_PerfectFrameVerticalStackedBar.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved exact-frame vertical runoff summary figure for {basin}\n")

print("All basin vertical plots successfully created with identical physical frames.")