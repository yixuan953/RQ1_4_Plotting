# This code is used to plot the crop production under different scenarios

# This code is used to plot the crop production under different scenarios
import os
import math
import pandas as pd
import matplotlib.pyplot as plt

# --- 1. Configuration ---
# Paths & Fractions (Adjust these paths if needed)
input_csv = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Impacts_NPCrop/Scenarios_Production_Runoff_Summary.csv"
current_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/1_Boundary" 
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig4"
os.makedirs(output_dir, exist_ok=True)

# N or P
# element = "N"
# dissolved_frac = 0.7  # Adjust if your actual fraction is less than 1.0
# var_name = "N [ktons]"

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

# Explicit scenario order (will be plotted from bottom to top)
scenario_order = [
    "Sus_Irri+Red_Fert",
    "Sus_Irri+Red_PFert", 
    "Sus_Irri+Red_NFert", 
    "Sus_Irri",
    "Baseline"   
]

crop_colors = {
    'Wheat': "#EBBB7C", 
    'Maize': "#B36767", 
    'Rice': "#C89CDC", 
    'Soybean': "#42949F"
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
        # We continue plotting the bars even if the boundary line file is missing

    # Pivot the data
    pivot_df = df_basin.pivot(index='Scenario', columns='Crop', values=var_name2)
    
    # Reindex to maintain logical sequence (Baseline at the bottom, Red_Fert at the top)
    pivot_df = pivot_df.reindex(scenario_order)
    pivot_df = pivot_df.reindex(columns=FinalCategories)

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
                color=crop_colors[crop],
                edgecolor='white', 
                linewidth=0.8,
                height=0.5,  
                alpha=0.9,
                zorder=3  # Puts bars above the grid lines
            )
            
            # Increment the left coordinate marker
            if left_x is None:
                left_x = lengths.fillna(0)
            else:
                left_x += lengths.fillna(0)

    # --- 5. Add Vertical Boundary Line ---
    if current_n_bound is not None:
        ax.axvline(
            x=current_n_bound, 
            color='black', 
            linestyle='-', 
            linewidth=2.0, 
            label='N Boundary',
            zorder=4  # Places line cleanly on top of the bars
        )
        print(f"  Added N boundary line at: {current_n_bound} ktons")

    # Clean Up Axes Visuals
    ax.spines['left'].set_color("#000000")
    ax.spines['bottom'].set_color("#000000")
    
    ax.set_yticklabels([])
    ax.tick_params(axis='x', labelsize=25)
    ax.grid(axis='x', linestyle='--', alpha=0.3, zorder=0)  # Grid lines on X-axis now
    ax.set_axisbelow(True)

    plt.tight_layout()
    
    # Save image
    out_path = os.path.join(output_dir, f"{basin}_Total_{element}_Runoff_HorizontalStackedBar.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved horizontal runoff summary figure for {basin}\n")

print("All basin horizontal plots successfully created.")