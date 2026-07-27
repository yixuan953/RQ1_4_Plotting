# This script is used to plot the required implementation to achieve the targeted footprint
# The plot contains 2 information
# 1 - What is the targeted footprint if a certain yield is needed (shadow area)
# 2 - How will footprint change with implementation of field management practices (or) and buffer zone (3 lines in total), range
# 3 - How will the targeted footrpint change if N, P delivery from point sources can be reduced to 0

import os
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- Configuration & Paths ---
basins = ["LaPlata", "Rhine", "Indus", "Yangtze"]
target_footprint_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Reconcile/NPFootprint_End-of-Pipe_Removal_Ratios.csv"

# Boundary Directory Paths
baseline_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/1_Boundary"
sewage50_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Reconcile/50Sewage"
sewage75_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Reconcile/75Sewage"
zero_point_source_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/3_TradeOffs"

output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig5/Required_Implementation_Rate.png"

# Yield scenarios for vertical gradient increments
yield_scenarios = [1, 0.9, 0.8, 0.7, 0.6, 0.5]

# X-axis scaling definitions
Implementation_Rates = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
tick_positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
tick_labels = ["0", "20", "40", "60", "80", "100"]

# --- Reduction Coefficients (with Ranges) ---
coefficients = {
    "N": {
        "FM": {"mid": 0.6, "min": 0.5, "max": 0.7},
        "BZ": {"mid": 0.5, "min": 0.2, "max": 0.8}
    },
    "P": {
        "FM": {"mid": 0.8, "min": 0.5, "max": 0.83},
        "BZ": {"mid": 0.67, "min": 0.1, "max": 1.1}
    }
}

# coefficients = {
#     "N": {
#         "FM": {"mid": 0.6, "min": 0.3, "max": 1.4},
#         "BZ": {"mid": 0.5, "min": 0.3, "max": 0.7}
#     },
#     "P": {
#         "FM": {"mid": 0.9, "min": 0.8, "max": 1.0},
#         "BZ": {"mid": 0.67, "min": 0.4, "max": 0.95}
#     }
# }

# --- Load Target Footprint Data ---
df = pd.read_csv(target_footprint_dir)

# --- Plot Setup ---
fig, axes = plt.subplots(2, 4, figsize=(24, 13), sharex="col", sharey=False)
elements = ["N", "P"]

shadow_cmap = mcolors.LinearSegmentedColormap.from_list(
    "yield_gradient", ["#ABCB87", "#6796CB"]
)

def get_nice_ticks(max_val):
    """Generates exactly 4 clean ticks (0, step, step*2, step*3) using standard intervals."""
    target_step = max_val / 3.1
    base_intervals = [0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    
    selected_step = base_intervals[-1]
    for interval in base_intervals:
        if interval >= target_step:
            selected_step = interval
            break
            
    ticks = [0.0, selected_step, selected_step * 2, selected_step * 3]
    ylim_max = selected_step * 3.3  
    return ticks, ylim_max

def get_cropland_boundary(folder, basin, element):
    file_path = os.path.join(folder, f"{basin}_boundaries_sum.csv")
    if not os.path.exists(file_path):
        return None
    try:
        b_df = pd.read_csv(file_path)
        label_col = b_df.columns[0]
        b_df = b_df.set_index(label_col)
        col_name = f"{element} [ktons]"
        return float(b_df.loc["Cropland", col_name])
    except Exception as e:
        print(f"Error reading boundary for {basin} {element} in {folder}: {e}")
        return None

for row_idx, elem in enumerate(elements):
    req_col = f"Required {elem} footprint"
    curr_col = f"{elem} footprint"
    coefs = coefficients[elem]
    
    for col_idx, basin in enumerate(basins):
        ax = axes[row_idx, col_idx]
        basin_max_val = 0.0

        # 1. Fetch all boundaries
        base_bound = get_cropland_boundary(baseline_boundary_dir, basin, elem)
        sewage50_bound = get_cropland_boundary(sewage50_boundary_dir, basin, elem)
        sewage75_bound = get_cropland_boundary(sewage75_boundary_dir, basin, elem)
        zero_bound = get_cropland_boundary(zero_point_source_boundary_dir, basin, elem)
        
        # Calculate scaling ratios relative to baseline boundary
        ratio_base = 1.0
        ratio_50 = (sewage50_bound / base_bound) if (base_bound and sewage50_bound) else 1.0
        ratio_75 = (sewage75_bound / base_bound) if (base_bound and sewage75_bound) else 1.0
        ratio_zero = (zero_bound / base_bound) if (base_bound and zero_bound) else 1.0

        # Arrays to hold values for the 4 distinct quarters
        zones_data = {
            "base": [],
            "50sewage": [],
            "75sewage": [],
            "zero": []
        }
        
        scenario_name = "Baseline yield + Required footprint"
        match = df[(df["Basin"] == basin) & (df["Scenario"] == scenario_name)]
        
        for yield_lvl in yield_scenarios:
            if not match.empty:
                val_base = match.iloc[0][req_col] * (1 / yield_lvl)
                
                v_base = val_base * ratio_base
                v_50 = val_base * ratio_50
                v_75 = val_base * ratio_75
                v_zero = val_base * ratio_zero
                
                zones_data["base"].append(v_base)
                zones_data["50sewage"].append(v_50)
                zones_data["75sewage"].append(v_75)
                zones_data["zero"].append(v_zero)
                
                basin_max_val = max(basin_max_val, v_base, v_50, v_75, v_zero)
            else:
                for key in zones_data:
                    zones_data[key].append(None)

        # 2. Extract Baseline Current Footprint
        baseline_scenario = "Baseline yield + Current footprint"
        baseline_match = df[(df["Basin"] == basin) & (df["Scenario"] == baseline_scenario)]
        current_fp = baseline_match.iloc[0][curr_col] if not baseline_match.empty else 0
        basin_max_val = max(basin_max_val, current_fp)

        # 3. Calculate Trajectories & Range Boundaries
        def calc_trajectory(c_val):
            return (current_fp * (1 - Implementation_Rates)) + (current_fp * Implementation_Rates * c_val)

        fp_with_FM = calc_trajectory(coefs["FM"]["mid"])
        fp_with_BZ = calc_trajectory(coefs["BZ"]["mid"])
        fp_with_Both = calc_trajectory(coefs["FM"]["mid"] * coefs["BZ"]["mid"])

        fp_with_FM_min = calc_trajectory(coefs["FM"]["min"])
        fp_with_FM_max = calc_trajectory(coefs["FM"]["max"])
        fp_with_BZ_min = calc_trajectory(coefs["BZ"]["min"])
        fp_with_BZ_max = calc_trajectory(coefs["BZ"]["max"])
        fp_with_Both_min = calc_trajectory(coefs["FM"]["min"] * coefs["BZ"]["min"])
        fp_with_Both_max = calc_trajectory(coefs["FM"]["max"] * coefs["BZ"]["max"])

        basin_max_val = max(basin_max_val, fp_with_FM_max.max(), fp_with_BZ_max.max(), fp_with_Both_max.max())

        # 4. Handle Subplot Framing Constraints
        y_ticks, y_limit_max = get_nice_ticks(basin_max_val)
        ax.set_ylim(0, y_limit_max)
        ax.set_yticks(y_ticks)
        
        y_tick_labels = [f"{t:g}" for t in y_ticks]
        ax.set_yticklabels(y_tick_labels, fontsize=28)

        # 5. Draw 4-Part Background Gradient Shadow Bands
        zone_bounds = [
            (0.00, 0.4, "base"),
            (0.4, 0.6, "50sewage"),
            (0.6, 0.8, "75sewage"),
            (0.8, 1.0, "zero")
        ]

        for x_start, x_end, key in zone_bounds:
            valid_vals = [(v, idx) for idx, v in enumerate(zones_data[key]) if v is not None]
            if len(valid_vals) > 1:
                # Highlight absolute bottom boundary (yield scenario 1)
                y_bottom = valid_vals[0][0]
                ax.plot([x_start, x_end], [y_bottom, y_bottom], color=shadow_cmap(0.0), linestyle="-", linewidth=3.5, zorder=2)
                
                # Highlight absolute top boundary (yield scenario 0.5)
                y_top = valid_vals[-1][0]
                ax.plot([x_start, x_end], [y_top, y_top], color=shadow_cmap(1.0), linestyle="-", linewidth=3.5, zorder=2)

                # Draw internal segments and transition grid lines
                for i in range(len(valid_vals) - 1):
                    y1, idx1 = valid_vals[i]
                    y2, _ = valid_vals[i+1]
                    color_frac = idx1 / float(len(yield_scenarios) - 1)
                    
                    ax.fill_between([x_start, x_end], min(y1, y2), max(y1, y2), 
                                    color=shadow_cmap(color_frac), alpha=0.32, zorder=0, linewidth=0)
                    ax.plot([x_start, x_end], [y1, y1], color="white", linestyle="-", linewidth=1.2, zorder=1)
                ax.plot([x_start, x_end], [valid_vals[-1][0], valid_vals[-1][0]], color="white", linestyle="-", linewidth=1.2, zorder=1)

        # Vertical sector dividers at 25%, 50%, and 75% positions
        for separator in [0.4, 0.6, 0.8, 1.0]:
            ax.axvline(x=separator, color="#95a5a6", linestyle=":", linewidth=1.5, alpha=0.6, zorder=2)

        # 6. Plot Main Implementation Lines & Cloud Ranges with Sharper/Obvious Edges
        # Field Management (FM) - High Z-Order values keep these completely on top
        ax.plot(Implementation_Rates, fp_with_FM, label="Field Management (FM)", 
                color="#FFA239", linestyle="--", linewidth=4.0, marker="s", markersize=11, zorder=6)
        ax.fill_between(Implementation_Rates, fp_with_FM_min, fp_with_FM_max, 
                        color="#B2B1B1", alpha=0.2, zorder=3)
        ax.plot(Implementation_Rates, fp_with_FM_min, color="#FFA239", linestyle="-", linewidth=2, alpha=0.6, zorder=4)
        ax.plot(Implementation_Rates, fp_with_FM_max, color="#FFA239", linestyle="-", linewidth=2, alpha=0.6, zorder=4)

        # Buffer Zone (BZ)
        ax.plot(Implementation_Rates, fp_with_BZ, label="Buffer Zone (BZ)", 
                color="#232F72", linestyle="--", linewidth=3.0, marker="o", markersize=9, zorder=7)
        ax.fill_between(Implementation_Rates, fp_with_BZ_min, fp_with_BZ_max, 
                        color="#B2B1B1", alpha=0.2, zorder=3)
        ax.plot(Implementation_Rates, fp_with_BZ_min, color="#232F72", linestyle="-", linewidth=2, alpha=0.6, zorder=4)
        ax.plot(Implementation_Rates, fp_with_BZ_max, color="#232F72", linestyle="-", linewidth=2, alpha=0.6, zorder=4)

        # Combined (FM + BZ)
        ax.plot(Implementation_Rates, fp_with_Both, label="FM + BZ Combined", 
                color="#C13383", linestyle="--", linewidth=3.5, marker="^", markersize=10, zorder=8)
        ax.fill_between(Implementation_Rates, fp_with_Both_min, fp_with_Both_max, 
                        color="#B2B1B1", alpha=0.2, zorder=3)
        ax.plot(Implementation_Rates, fp_with_Both_min, color="#C13383", linestyle="-", linewidth=2, alpha=0.6, zorder=4)
        ax.plot(Implementation_Rates, fp_with_Both_max, color="#C13383", linestyle="-", linewidth=2, alpha=0.6, zorder=4)

        # Baseline point anchor
        ax.plot(0.0, current_fp, marker="o", color="#4b5151", markersize=14, zorder=9)

        # --- Axis Polish ---
        ax.set_xlim(-0.06, 1.06)  # Retains frame margin padding
        ax.grid(True, linestyle=":", alpha=0.25, zorder=2)
        
        if row_idx == 1:
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, fontsize=30)
        else:
            ax.set_xticklabels([])

plt.tight_layout(rect=[0, 0.04, 1, 0.94])

# Save Output Configuration
os.makedirs(os.path.dirname(output_dir), exist_ok=True)
plt.savefig(output_dir, dpi=300)
print(f"Plot successfully saved to: {output_dir}")