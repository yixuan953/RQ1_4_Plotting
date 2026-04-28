import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mtick

# --- Configuration ---
Studyareas = ["Indus", "LaPlata", "Yangtze", "Rhine"]
# Mapping internal names to display names
Crop_Mapping = {
    "winterwheat": "Wheat",
    "maize": "Maize",
    "mainrice": "Rice",
    "soybean": "Soybean"
}
Categories = ["N & P", "P", "N"]
colors = ["#FF875F", "#FDD804", "#E41F7E"]  # Gray, Orange, Blue matching your excel chart

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Statistics/3_Exceedance"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig4"
os.makedirs(output_dir, exist_ok=True)

for basin in Studyareas:
    file_path = os.path.join(input_dir, f"{basin}_boundaries_ExceedanceHA.csv")
    if not os.path.exists(file_path):
        print(f"File not found for {basin}, skipping...")
        continue

    # Load data
    df = pd.read_csv(file_path)
    
    # Filter for the specific crops and maintain order
    df = df[df['Crop type'].isin(Crop_Mapping.keys())].copy()
    df['Display Name'] = df['Crop type'].map(Crop_Mapping)
    
    # Reverse order so Wheat is at the top in a horizontal plot
    df = df.iloc[::-1]
    
    num_crops = len(df)
    if num_crops == 0:
        continue

    # --- Plotting ---
    # Adjust figure height dynamically: ~2 inches per crop + padding
    fig_height = max(4, num_crops * 1.5)
    fig, ax = plt.subplots(figsize=(5, fig_height))

    y = np.arange(num_crops)
    width = 0.3  # width of the bars

    # Plot bars for N, P, and N & P
    # We multiply by 100 to convert decimal (e.g., 0.51) to percentage (51%)
    ax.barh(y - width, df['N & P'] * 100, width, label='N & P', color=colors[0], alpha = 0.8, edgecolor='black')
    ax.barh(y, df['P'] * 100, width, label='P', color=colors[1], alpha = 0.8, edgecolor='black')
    ax.barh(y + width, df['N'] * 100, width, label='N', color=colors[2], alpha = 0.8, edgecolor='black')

    # Formatting
    ax.xaxis.tick_top() # Moves the ticks/labels to the top
    ax.xaxis.set_label_position('top') # Moves the "Percentage..." title to the top
    ax.set_yticks(y)
    ax.set_yticklabels(df['Display Name'], fontsize=25)
    
    ax.set_xlim(0, 100)
    ax.tick_params(axis='both', labelsize=25)

    # Add grid lines for readability
    ax.grid(axis='x', linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    # Optional: Hide the bottom spine since the axis is on top
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    
    save_path = os.path.join(output_dir, f"{basin}_Exceedance_HA_bar.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved exceedance plot for {basin}")