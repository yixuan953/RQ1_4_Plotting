# This code is used to plot the crop production under different scenarios

import os
import pandas as pd
import matplotlib.pyplot as plt

# --- 1. Configuration ---
Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
FinalCategories = ["Soybean", "Rice", "Maize", "Wheat"]

# Explicit scenario order (will be plotted from bottom to top)
scenario_order = [
    "Water+NP",
    "Water+P",
    "Water+N",
    "Water",
    "Current"   
]

crop_hatches = {
    'Wheat': '|',     # Single forward slash (diagonal lines)
    'Maize': 'x',    # Single backslash (opposite diagonal lines)
    'Rice': '.',      # Single dot (much cleaner than '...')
    'Soybean': '\\'    # Single plus (simple grid, spacious)
}

scenario_face_colors = {
    "Water+NP": "#02B585",  # green
    "Water+P" : "#E8A921",  # orange
    "Water+N" : "#D970AA",  # purple-pink
    "Water"   : "#56B4E9",  # sky blue
    "Current" : "#999999"   # gray
}

input_csv = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/2_NPYield_Change/RespectBoundaries_Yield_Runoff_Summary.csv"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig3"
os.makedirs(output_dir, exist_ok=True)

# --- 2. Load and Clean Data ---
df = pd.read_csv(input_csv)

# Filter for only target crops
df = df[df['Crop'].isin(FinalCategories)]

# --- 3. Loop and Plot for Each Basin ---
for basin in Studyareas:
    df_basin = df[df['Basin'] == basin].copy()
    
    if df_basin.empty:
        print(f"No data found for basin: {basin}, skipping...")
        continue
        
    # Pivot the data
    pivot_df = df_basin.pivot(index='Scenario', columns='Crop', values='Production_ktons')
    
    # Reindex to maintain logical sequence (Baseline at the bottom, Red_Fert at the top)
    pivot_df = pivot_df.reindex(scenario_order)
    pivot_df = pivot_df.reindex(columns=FinalCategories)
    pivot_df = pivot_df / 1000.0

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Track baseline left position for horizontal stacking
    left_x = None 
    
    for crop in FinalCategories:
        if crop in pivot_df.columns:
            lengths = pivot_df[crop]
            
            # Draw horizontal stacked components using barh
            ax.barh(
                pivot_df.index,
                lengths,
                left=left_x,
                label=crop,
                hatch=crop_hatches[crop],
                facecolor=[scenario_face_colors.get(s, 'lightgray')
                            for s in pivot_df.index],
                edgecolor='black',
                linewidth=0.8,
                height=0.5,
                alpha=0.9
            )
                        
            # Increment the left coordinate marker
            if left_x is None:
                left_x = lengths.fillna(0)
            else:
                left_x += lengths.fillna(0)

    # Clean Up Axes Visuals
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color("#000000")
    ax.spines['bottom'].set_color("#000000")
    
    ax.set_yticklabels([])
    ax.tick_params(axis='x', labelsize=25)
    ax.grid(axis='x', linestyle='--', alpha=0.3, zorder=0)  # Grid lines on X-axis now
    ax.set_axisbelow(True)

    plt.tight_layout()
    
    # Save image
    out_path = os.path.join(output_dir, f"{basin}_Total_Production_HorizontalStackedBar.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved horizontal production summary figure for {basin}")

print("All basin horizontal plots successfully created.")