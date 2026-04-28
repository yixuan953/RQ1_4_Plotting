# This script is used to compare the N and P runoff under each scenarios with the boundaries

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import FixedLocator, FixedFormatter

# --- Paths ---
input_runoff_file = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/4_Red_Prod_Runoff/Production_Runoff_Summary.csv"
input_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/1_Boundary_load"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig4b_Runoff_Boundary"

os.makedirs(output_dir, exist_ok=True)

STUDY_AREAS = ["LaPlata", "Rhine", "Indus", "Yangtze"]
SCENARIOS = ['Baseline', 'Sus_Irri', 'Sus_Irri+Red_Fert'] 
CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}

def format_reduction_label(abs_red, pct_red):
    if pct_red < 0.1:
        return "(0%)"
    pct_str = f"{pct_red:.1f}%" if pct_red < 1.0 else f"{pct_red:.0f}%"
    abs_str = f"{abs_red:.2f} kt" if abs_red < 0.1 else f"{abs_red:.1f} kt"
    return f"-{abs_str} (-{pct_str})" # Added newline for better fit

def generate_butterfly_runoff():
    df_runoff = pd.read_csv(input_runoff_file)
    
    # Increased figsize for vertical clarity
    fig, (ax_n, ax_p) = plt.subplots(1, 2, figsize=(12, 40), sharey=True)
    # Balanced wspace to allow text without overlapping the middle
    plt.subplots_adjust(wspace=0.35)

    bar_height = 0.6
    scen_colors = ["#F5E9D8", '#2FA4D7', '#E76F2E'] 
    
    current_y = 0
    ytick_pos, ytick_labels = [], []

    for basin in STUDY_AREAS:
        boundary_path = os.path.join(input_boundary_dir, f"{basin}_crit_load_sum_updated.csv")
        if not os.path.exists(boundary_path): continue
        df_bound = pd.read_csv(boundary_path, index_col=0)
        basin_data = df_runoff[df_runoff['Basin'] == basin]
        
        for crop_group, subcrops in CROP_MAP.items():
            available_bounds = df_bound.index.intersection(subcrops)
            if available_bounds.empty: continue

            base_row = basin_data[(basin_data['Crop'] == crop_group) & (basin_data['Scenario'] == 'Baseline')]
            if base_row.empty: continue
            
            n_base_val = base_row['N_Runoff_ktons'].values[0]
            p_base_val = base_row['P_Runoff_ktons'].values[0]
            n_lim = df_bound.loc[available_bounds, 'N [ktons]'].sum()
            p_lim = df_bound.loc[available_bounds, 'P [ktons]'].sum()

            ytick_pos.append(current_y + (len(SCENARIOS)*bar_height)/2)
            ytick_labels.append(crop_group)

            for i, scen in enumerate(SCENARIOS):
                row = basin_data[(basin_data['Crop'] == crop_group) & (basin_data['Scenario'] == scen)]
                if row.empty: continue
                scen_row = row.iloc[0]
                y_pos = current_y + (i * bar_height)
                
                # We use 0.1 as the minimum for the log scale bars
                n_rat = max(scen_row['N_Runoff_ktons'] / n_lim, 0.1) if n_lim > 0 else 0.1
                p_rat = max(scen_row['P_Runoff_ktons'] / (p_lim * 0.25), 0.1) if p_lim > 0 else 0.1

                ax_n.barh(y_pos, n_rat, bar_height, color=scen_colors[i], edgecolor='white')
                ax_p.barh(y_pos, p_rat, bar_height, color=scen_colors[i], edgecolor='white')

                if scen != 'Baseline':
                    # N Logic: Place text further LEFT of the bar end
                    n_abs, n_pct = n_base_val - scen_row['N_Runoff_ktons'], ((n_base_val - scen_row['N_Runoff_ktons'])/n_base_val*100)
                    ax_n.text(n_rat / 1.15, y_pos, format_reduction_label(n_abs, n_pct), 
                              va='center', ha='right', fontsize=25, color="#34495e")
                    
                    # P Logic: Place text further RIGHT of the bar end
                    p_abs, p_pct = p_base_val - scen_row['P_Runoff_ktons'], ((p_base_val - scen_row['P_Runoff_ktons'])/p_base_val*100)
                    ax_p.text(p_rat * 1.15, y_pos, format_reduction_label(p_abs, p_pct), va='center', ha='left', fontsize=25, color='#34495e')

            current_y += (len(SCENARIOS) * bar_height) + 1.2 # Increased crop gap
        current_y += 2.0 # Increased basin gap

    # Log scale ticks
    custom_ticks = [0.1, 1, 2, 5, 10, 50]
    custom_labels = ["0", "1", "2", "5", "10", "50"]

    for ax, title in zip([ax_n, ax_p], ['N Runoff / Limit', 'P Runoff / Limit']):
        ax.set_xscale('log')
        ax.set_xlim(0.1, 70) # Increased limit to 70 to give text breathing room
        ax.xaxis.set_major_locator(FixedLocator(custom_ticks))
        ax.xaxis.set_major_formatter(FixedFormatter(custom_labels))
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(title, fontsize=25, labelpad=20)
        ax.tick_params(axis='x', labelsize=25)
        
        # Frame Logic
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(True)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.axvline(1.0, color='red', linestyle='--', linewidth=2, alpha=0.4, zorder=0)

    # N-specific center spine
    ax_n.spines['right'].set_visible(True)
    ax_n.spines['right'].set_linewidth(1.5)
    ax_n.invert_xaxis()
    ax_n.yaxis.set_ticks_position('none')
    
    # P-specific center spine
    ax_p.spines['left'].set_visible(True)
    ax_p.spines['left'].set_linewidth(1.5)

    # Centering Crop Labels in the gutter
    ax_n.set_yticks(ytick_pos)
    ax_n.set_yticklabels([]) 

    plt.gca().invert_yaxis()
    plt.savefig(os.path.join(output_dir, "Butterfly_Runoff_Final_Improved.png"), dpi=300, bbox_inches='tight')
    plt.show()

generate_butterfly_runoff()