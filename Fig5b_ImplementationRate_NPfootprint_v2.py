# This script is used to plot the required implementation to achieve the targeted footprint
# The plot contains 2 information
# 1 - What is the targeted footprint if a certain yield is needed (shadow area)
# 2 - How will footprint change with implementation of field management practices (or) and buffer zone (3 lines in total)


import os
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- Configuration & Paths ---
basins = ["LaPlata", "Rhine", "Indus", "Yangtze"]
target_footprint_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Reconcile/NPFootprint_End-of-Pipe_Removal_Ratios.csv"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig5/Required_Implementation_Rate.png"

# The exact yield scenarios to extract for the background gradient limits
yield_scenarios = [1, 0.9, 0.8, 0.7, 0.6, 0.5]

# X-axis scale definitions for implementation rates
Implementation_Rates = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
tick_positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
tick_labels = ["0", "20", "40", "60", "80", "100"]

# --- Reduction Coefficients ---
N_FootPrint_after_FieldManagement = 0.4  # 40% of current
P_FootPrint_after_FieldManagement = 0.8  # 80% of current

N_FootPrint_after_BufferZone = 0.4  # 40% of current
P_FootPrint_after_BufferZone = 0.7  # 70% of current

# --- Load Data ---
df = pd.read_csv(target_footprint_dir)

# --- Plot Setup ---
fig, axes = plt.subplots(2, 4, figsize=(24, 11), sharex="col", sharey="row")
elements = ["N", "P"]

# Custom continuous colormap running precisely from requested Green (#84994F) to Red (#B45253)
shadow_cmap = mcolors.LinearSegmentedColormap.from_list(
    "yield_gradient", ["#9FC04A", "#D04648"]
)

for row_idx, elem in enumerate(elements):
    req_col = f"Required {elem} footprint"
    curr_col = f"{elem} footprint"

    # Assign management coefficients per element
    fp_after_FM = N_FootPrint_after_FieldManagement if elem == "N" else P_FootPrint_after_FieldManagement
    fp_after_BZ = N_FootPrint_after_BufferZone if elem == "N" else P_FootPrint_after_BufferZone
    fp_after_Both = fp_after_FM * fp_after_BZ

    # Track maximum value to adjust Y-limits dynamically per row
    row_max_val = 0.0

    for col_idx, basin in enumerate(basins):
        ax = axes[row_idx, col_idx]

        # 1. Extract the Targeted Footprint values for each yield scenario
        target_fps = []
        for yield_lvl in yield_scenarios:
            scenario_name = f"Baseline yield + Required footprint"
            match = df[(df["Basin"] == basin) & (df["Scenario"] == scenario_name)]
            
            if not match.empty:
                val = match.iloc[0][req_col] * (1/yield_lvl)
                target_fps.append(val)
                row_max_val = max(row_max_val, val)
            else:
                target_fps.append(None)

        # 2. Draw the changing color shadow bands based on the REAL target footprint positions
        # Filter out any missing/None values while keeping track of their original index for color accuracy
        valid_fps = [(val, idx) for idx, val in enumerate(target_fps) if val is not None]
        
        if len(valid_fps) > 1:
            for i in range(len(valid_fps) - 1):
                y1, idx1 = valid_fps[i]
                y2, idx2 = valid_fps[i+1]
                
                ymin_band = min(y1, y2)
                ymax_band = max(y1, y2)
                
                # Color fraction is dynamically tied to the actual yield reduction scenario index
                # (Baseline = 0.0 -> Green, 50% Reduction = 1.0 -> Red)
                color_frac = idx1 / float(len(yield_scenarios) - 1)
                facecolor = shadow_cmap(color_frac)
                
                # ax.axhspan completely fills from left to right edge using the real data values
                ax.axhspan(
                    ymin=ymin_band,
                    ymax=ymax_band,
                    color=facecolor,
                    alpha=0.35, # Soft transparency for a clean, professional aesthetic
                    zorder=0
                )
                
                # Draw sharp white divider lines at each real footprint level boundary
                ax.axhline(y=y1, color="white", linestyle="-", linewidth=1.0, zorder=1)
                
            # Draw the final top boundary line for the last scenario row
            ax.axhline(y=valid_fps[-1][0], color="white", linestyle="-", linewidth=1.0, zorder=1)

        # 3. Extract Baseline Current Footprint to calculate management curves
        baseline_scenario = "Baseline yield + Current footprint"
        baseline_match = df[(df["Basin"] == basin) & (df["Scenario"] == baseline_scenario)]
        current_fp = baseline_match.iloc[0][curr_col] if not baseline_match.empty else 0
        row_max_val = max(row_max_val, current_fp)

        # 4. Calculate Management Curve Trajectories over implementation rates
        fp_with_FM = (current_fp * (1 - Implementation_Rates)) + (current_fp * Implementation_Rates * fp_after_FM)
        fp_with_BZ = (current_fp * (1 - Implementation_Rates)) + (current_fp * Implementation_Rates * fp_after_BZ)
        fp_with_Both = (current_fp * (1 - Implementation_Rates)) + (current_fp * Implementation_Rates * fp_after_Both)

        # 5. Plot Implementation Curves

        # Scenario trajectories plotted explicitly from 10% to 100% implementation rate
        ax.plot(
            Implementation_Rates[0:],
            fp_with_FM[0:],
            label="Field Management (FM)",
            color="#34495e",
            linestyle="--",
            linewidth=4.5,
            marker="s",
            markersize=11,
            zorder=3,
        )
        ax.plot(
            Implementation_Rates[0:],
            fp_with_BZ[0:],
            label="Buffer Zone (BZ)",
            color="#16a98c",
            linestyle="--",
            linewidth=3.3,
            marker="o",
            markersize=8,
            zorder=4,
        )
        ax.plot(
            Implementation_Rates[0:],
            fp_with_Both[0:],
            label="FM + BZ Combined",
            color="#c78a2d",
            linestyle="--",
            linewidth=4,
            marker="^",
            markersize=10,
            zorder=5,
        )

        # Single grey dot representing the initial current footprint state
        ax.plot(
            0.0, 
            current_fp, 
            marker="o", 
            color="#5D6161", 
            markersize=10, 
            zorder=6, 
            label="Current Footprint"
        )
        

        # --- Subplot Formatting ---
        # Expanded X-limits slightly so boundary dots are not clipped at the edges
        ax.set_xlim(-0.08, 1.08)
        ax.grid(True, linestyle=":", alpha=0.3, zorder=2)
        
        # Clean up tick management to avoid internal label clustering
        if row_idx == 1:
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, fontsize=30)
        else:
            ax.set_xticklabels([])

    # Apply row-wide uniform Y-limits with 15% breathing room
    for col_idx in range(4):
        axes[row_idx, col_idx].set_ylim(0, row_max_val * 1.15)
        # Use MaxNLocator to explicitly reduce the density of vertical ticks
        axes[row_idx, col_idx].yaxis.set_major_locator(plt.MaxNLocator(nbins=4))
        # Enhanced font size configuration for y-axis tick marks
        axes[row_idx, col_idx].tick_params(axis='y', labelsize=30)


plt.tight_layout(rect=[0, 0.05, 1, 0.92])

# Save Plot
os.makedirs(os.path.dirname(output_dir), exist_ok=True)
plt.savefig(output_dir, dpi=300)
print(f"Plot successfully saved to: {output_dir}")