import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import PercentFormatter

# Define paths
input_file = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Reconcile/NPFootprint_End-of-Pipe_Removal_Ratios.csv"
output_image = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig5/Removal_Ratio_TradeOffs.png"

if not os.path.exists(input_file):
    raise FileNotFoundError(f"Please run the updated calculation script first to generate: {input_file}")

df = pd.read_csv(input_file)

# Parse scenario string metadata
def parse_scenario(row):
    parts = row["Scenario"].split(" yield + ")
    yield_scen = parts[0]
    footprint_part = parts[1]
    
    if "Current footprint" in footprint_part:
        return yield_scen, "0%"
    elif "Required footprint" in footprint_part:
        return yield_scen, "Required"
    else:
        pct = footprint_part.split(" ")[0]
        return yield_scen, pct

df[["Yield_Scenario", "Footprint_Reduction"]] = df.apply(parse_scenario, axis=1, result_type="expand")

# Filter out the "Required" footprints to isolate the clean 0% to 80% reduction line trend
df_plot = df[df["Footprint_Reduction"] != "Required"].copy()

# Map percentage strings to structural numerical x-axis values
pct_map = {f"{i}%": i for i in [0, 10, 20, 30, 40, 50, 60, 70, 80]}
df_plot["X_val"] = df_plot["Footprint_Reduction"].map(pct_map)
df_plot = df_plot.sort_values(by="X_val")

# Structural grid layout configurations
basins = ["LaPlata", "Rhine", "Indus", "Yangtze"]
yield_order = ["Baseline", "Sus_Irrigation", "5% Reduction", "10% Reduction", "20% Reduction", "30% Reduction", "40% Reduction", "50% Reduction"]

# Configure custom colors and transparency levels (alphas)
base_color_green = "#26A382"
color_dict = {
    "Baseline": "grey",
    "Sus_Irrigation": "#4793BE",
    "5% Reduction": base_color_green,
    "10% Reduction": base_color_green,
    "20% Reduction": base_color_green,
    "30% Reduction": base_color_green,
    "40% Reduction": base_color_green,
    "50% Reduction": base_color_green
}

alpha_dict = {
    "Baseline": 1.0,
    "Sus_Irrigation": 1.0,
    "5% Reduction": 0.3,
    "10% Reduction": 0.4,
    "20% Reduction": 0.5,
    "30% Reduction": 0.6,
    "40% Reduction": 0.8,
    "50% Reduction": 1.0
}

# Initialize figure layout with larger default text elements
plt.rcParams.update({'font.size': 14}) 
fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 22), sharex=True, constrained_layout=True)

for row_idx, basin in enumerate(basins):
    df_basin = df_plot[df_plot["Basin"] == basin]
    
    ax_n = axes[row_idx, 0]
    ax_p = axes[row_idx, 1]
    
    for yield_scen in yield_order:
        df_sub = df_basin[df_basin["Yield_Scenario"] == yield_scen]
        if df_sub.empty:
            continue
            
        color = color_dict[yield_scen]
        alpha = alpha_dict[yield_scen]
        
        # Plot N dynamic lines
        ax_n.plot(df_sub["X_val"], df_sub["N Boundary removal ratio"], 
                  marker='o', linestyle='-', linewidth=3, markersize=8,
                  color=color, alpha=alpha, label=yield_scen)
        
        # Plot P dynamic lines
        ax_p.plot(df_sub["X_val"], df_sub["P Boundary removal ratio"], 
                  marker='s', linestyle='-', linewidth=3, markersize=8,
                  color=color, alpha=alpha, label=yield_scen)
        
    # # Formatting subplots with scaled-up font labels
    # ax_n.set_ylabel(f"{basin}\nN Removal Ratio", fontsize=16, fontweight='bold')
    # ax_p.set_ylabel(f"P Removal Ratio", fontsize=16, fontweight='bold')
    
    ax_n.grid(True, linestyle=":", alpha=0.6)
    ax_p.grid(True, linestyle=":", alpha=0.6)
    
    ax_n.set_ylim(-0.05, 1.05)
    ax_p.set_ylim(-0.05, 1.05)
    
    # Format Left Y-axis (N) to display cleanly as percentages
    ax_n.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    
    # Format Right Y-axis (P) to hide numeric text labels while maintaining layout sizing
    ax_p.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax_n.tick_params(axis='both', labelsize=25)
    ax_p.tick_params(axis='both', labelsize=25)

# # Header Titles Styling
# axes[0, 0].set_title("Nitrogen (N) Removal Dynamics", fontsize=18, fontweight='bold', pad=15)
# axes[0, 1].set_title("Phosphorus (P) Removal Dynamics", fontsize=18, fontweight='bold', pad=15)

x_ticks = [0, 20, 40, 60, 80]
x_labels = ["0%", "20%", "40%", "60%", "80%"]

for col_idx in [0, 1]:
    axes[3, col_idx].set_xticks(x_ticks)
    axes[3, col_idx].set_xticklabels(x_labels, rotation=30, ha='right', fontsize=25)
    # axes[3, col_idx].set_xlabel("Footprint Reduction (%)", fontsize=16, fontweight='bold', labelpad=12)

# # Unified global legend located safely underneath subplots
# handles, labels = axes[0, 0].get_legend_handles_labels()
# fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.04), 
#            ncol=4, fontsize=14, title="Yield Scenarios", title_fontsize=16, frameon=True)

# Export image output
plt.savefig(output_image, dpi=300, bbox_inches='tight')
plt.close()

print(f"Figure updated successfully. Layout preserved with hidden P labels and percentage scales. Saved to: {output_image}")