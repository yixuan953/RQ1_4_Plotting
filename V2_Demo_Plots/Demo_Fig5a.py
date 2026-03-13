import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D

# --- Configuration ---
basins = ["LaPlata", "Rhine", "Indus", "Yangtze"]
crops = ["Soybean", "Wheat", "Rice", "Maize"]
systems = ["Irrigated", "Rainfed"]
highlight_scenarios = ["Baseline", "Red_prop"]
highlight_colors = {"Baseline": "orange", "Red_prop": "#1b9e77"}

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/4_Fertilization_Red/4_2_Summary"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig5/5a"
os.makedirs(output_dir, exist_ok=True)

# Set global font sizes for high visibility
plt.rcParams.update({'font.size': 26, 'axes.labelsize': 26, 'axes.titlesize': 26, 
                     'xtick.labelsize': 26, 'ytick.labelsize': 26, 'legend.fontsize': 24})
sns.set_style("ticks")

crop_colors = {"Wheat": "#B379B4", "Rice": "#4D97CD", "Maize": "#eb8b39", "Soybean": "#72b285"}
system_markers = {"Irrigated": "o", "Rainfed": "s"}

def masked_stats(x, stat_type):
    clean_x = x[x != 0].dropna()
    if clean_x.empty: return 0.0
    if stat_type == 'mean': return clean_x.mean()
    if stat_type == 'p25':  return np.percentile(clean_x, 25)
    if stat_type == 'p75':  return np.percentile(clean_x, 75)
    return 0.0

for basin in basins:
    avg_file = os.path.join(input_dir, f"{basin}_basin_Avg_Summary.csv")
    total_file = os.path.join(input_dir, f"{basin}_basin_Total_Summary.csv")
    if not (os.path.exists(avg_file) and os.path.exists(total_file)): continue

    df_avg_raw = pd.read_csv(avg_file)
    df_total_raw = pd.read_csv(total_file)

    avg_summary = df_avg_raw.groupby(['Scenario', 'Crop', 'System']).agg({'Fert': lambda x: masked_stats(x, 'mean')}).reset_index()
    total_summary = df_total_raw.groupby(['Scenario', 'Crop', 'System']).agg({
        'N_runoff': [lambda x: masked_stats(x, 'mean'), lambda x: masked_stats(x, 'p25'), lambda x: masked_stats(x, 'p75')],
        'Yield': [lambda x: masked_stats(x, 'mean'), lambda x: masked_stats(x, 'p25'), lambda x: masked_stats(x, 'p75')]
    }).reset_index()

    total_summary.columns = ['Scenario', 'Crop', 'System', 'N_runoff_mean', 'N_runoff_p25', 'N_runoff_p75', 'Yield_mean', 'Yield_p25', 'Yield_p75']
    df_plot = pd.merge(avg_summary, total_summary, on=['Scenario', 'Crop', 'System'])
    
    # Unit conversion and Ratio calculation
    cols_to_scale = ['N_runoff_mean', 'N_runoff_p25', 'N_runoff_p75', 'Yield_mean', 'Yield_p25', 'Yield_p75']
    df_plot[cols_to_scale] = df_plot[cols_to_scale] * 0.000001 
    df_plot['Yield_N_Ratio_mean'] = np.where(df_plot['N_runoff_mean'] > 0, 
                                            df_plot['Yield_mean'] / df_plot['N_runoff_mean'], 0)
    # p25 Ratio (Lower efficiency: Low Yield / High Runoff)
    df_plot['Yield_N_Ratio_p25'] = np.where(df_plot['N_runoff_p75'] > 0, 
                                           df_plot['Yield_p25'] / df_plot['N_runoff_p75'], 0)
    # p75 Ratio (Higher efficiency: High Yield / Low Runoff)
    df_plot['Yield_N_Ratio_p75'] = np.where(df_plot['N_runoff_p25'] > 0, 
                                           df_plot['Yield_p75'] / df_plot['N_runoff_p25'], 0)

    # 1. Filter out all "Inc_" scenarios
    df_plot = df_plot[~df_plot['Scenario'].str.startswith('Inc_')].copy()
    threshold_lookup = {
            "Yangtze": 0.75,
            "Indus": 0.75,
            "Rhine": 0.65,
            "LaPlata": 0.20
        }
    for crop in crops:
        df_crop = df_plot[df_plot['Crop'] == crop].copy()
        if df_crop.empty: continue

        fig, axes = plt.subplots(3, 1, figsize=(12, 22), sharex=True)

        metrics = [
                    ('N_runoff_mean', 'N_runoff_p25', 'N_runoff_p75', 'Total N Runoff (ktons N)'),
                    ('Yield_mean', 'Yield_p25', 'Yield_p75', 'Total Yield (ktons)'),
                    ('Yield_N_Ratio_mean', 'Yield_N_Ratio_p25', 'Yield_N_Ratio_p75', 'Yield / N Runoff Ratio')
                ]

        for i, (mean_col, p25_col, p75_col, ylabel) in enumerate(metrics):
            ax = axes[i]
            # --- ADD OR REPLACE HERE ---
            ax.set_ylabel(ylabel, labelpad=15)
            # Customized grey gridlines
            ax.grid(True, which='both', axis='y', linestyle='--', color='grey', alpha=0.2)
            ax.grid(True, which='major', axis='x', linestyle='-', color='grey', alpha=0.2)

            for system in systems:
                # 4. Custom System Filtering logic
                if basin == "LaPlata" and system == "Irrigated": continue
                if basin == "Rhine" and crop == "Wheat" and system == "Irrigated": continue
                if basin == "Yangtze" and crop == "Rice" and system == "Rainfed": continue

                subset = df_crop[df_crop['System'] == system].copy()
                if subset.empty: continue
                subset = subset.sort_values('Fert')
                color = crop_colors.get(crop, 'blue')

                # --------- Plot Threshold Line on the Second Subplot (Yield) ---
                if i == 1: # Index 1 is the Yield subplot
                    baseline_yield = subset[subset['Scenario'] == 'Baseline']['Yield_mean'].values
                    if len(baseline_yield) > 0:
                        pct = threshold_lookup.get(basin, 1.0)
                        threshold_val = baseline_yield[0] * pct
                        
                        # Draw horizontal line
                        ax.axhline(y=threshold_val, color="red", 
                                   linestyle=':', linewidth=3, alpha=0.8, zorder=1)
                        
                # 2. Range Plotting: Only Baseline and Red_xx (Exclude Red_prop)
                if p25_col and p75_col:
                    subset_range = subset[(subset['Scenario'] == 'Baseline') | 
                                          (subset['Scenario'].str.startswith('Red_') & (subset['Scenario'] != 'Red_prop'))].copy()
                    ax.fill_between(subset_range['Fert'], subset_range[p25_col], subset_range[p75_col], color=color, alpha=0.15, lw=0)
                
                # Dashed Line Connectivity (Exclude Red_prop)
                subset_line = subset[subset['Scenario'] != 'Red_prop'].copy()
                ax.plot(subset_line['Fert'], subset_line[mean_col], color=color, linestyle='--', lw=2.5, alpha=0.7)
                
                # Scatter background dots
                ax.scatter(subset['Fert'], subset[mean_col], color=color, marker=system_markers[system], s=100, alpha=0.4)

                # 3. Highlight Baseline and Red_prop with "Boxplot" style
                for sc in highlight_scenarios:
                    h_point = subset[subset['Scenario'] == sc]
                    if not h_point.empty:
                        # Circle/Square + Edge line highlight
                        ax.scatter(h_point['Fert'], h_point[mean_col], 
                                   color=highlight_colors[sc], marker=system_markers[system], s=350, 
                                   edgecolor='white', linewidth=1, zorder=12)

            ax.set_ylabel(ylabel, labelpad=15)
            ax.grid(True, linestyle=':', alpha=0.4)

        axes[-1].set_xlabel("Average N fertilizer input (kg N/ha)", labelpad=10)

        # Legend
        system_legends = [Line2D([0], [0], marker=system_markers[s], color='gray', label=s, linestyle='None', markersize=12) for s in systems]
        scen_legends = [Line2D([0], [0], marker='s', color=highlight_colors[sc], markeredgecolor='black', 
                               markersize=15, label=sc, linestyle='None') for sc in highlight_scenarios]
        
        fig.legend(handles=system_legends + scen_legends, loc='lower center', bbox_to_anchor=(0.5, 0.02), ncol=4, frameon=True)
        plt.tight_layout(rect=[0, 0.07, 1, 0.95])
        
        plt.savefig(os.path.join(output_dir, f"{basin}_{crop}_Summary_Vertical.png"), dpi=300)
        plt.close()