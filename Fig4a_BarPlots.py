# This script is used to plot the bar plots for changes in:
# 1 - Crop production
# 2 - Cropland N runoff
# 3 - Cropland P runoff

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# --- Configuration ---
Boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/1_Boundary_load"
Prod_Runoff_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/4_Red_Prod_Runoff"
figure_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig4_Crop_Prod_Red"
os.makedirs(figure_dir, exist_ok=True)

STUDY_AREAS = ["Indus", "LaPlata", "Yangtze", "Rhine"]
all_crops = ['Wheat', 'Maize', 'Rice', 'Soybean']
SCENARIOS = ["Baseline", "Sus_Irri", "Sus_Irri+Red_Fert"]
COLORS = ["#7f8c8d", "#5bb5f0", "#f8ad6b", "#a4f3c5"] # Baseline, Sus_Irr, Red_Fert, Inc_Fert
PP_Prop = 0.25

# --- Helper Functions ---

def Get_Boundaries(basin, crop_display_name):
    file_path = os.path.join(Boundary_dir, f"{basin}_crit_load_sum_updated.csv")
    if not os.path.exists(file_path):
        print(f"WARNING: Boundary file not found for {basin}")
        return None, None
    df_b = pd.read_csv(file_path, index_col=0)

    mapping = {
        "Wheat": ["winterwheat"],
        "Maize": ["maize"],
        "Rice": ["mainrice", "secondrice"],
        "Soybean": ["soybean"]
    }
    
    target_keys = mapping.get(crop_display_name, [])
    n_sum = 0
    p_sum = 0
    found_any = False

    for key in target_keys:
        if key in df_b.index:
            n_sum += df_b.loc[key, "N [ktons]"]
            p_sum += df_b.loc[key, "P [ktons]"]
            found_any = True
            
    if not found_any:
        return None, None
    p_sum_final = p_sum * PP_Prop
    
    return n_sum, p_sum_final

def get_smart_ticks(max_val):
    """Determines a nice rounding limit for the Y-axis."""
    if max_val <= 1: return 1, 0.2
    if max_val <= 5: return 5, 1
    if max_val <= 10: return 10, 2
    if max_val <= 20: return 20, 5
    if max_val <= 50: return 50, 10
    if max_val <= 100: return 100, 20
    if max_val <= 200: return 200, 50
    if max_val <= 500: return 500, 100
    if max_val <= 1000: return 1000, 200
    upper = np.ceil(max_val / 100) * 100
    return upper, upper / 5

def Plot_Basin_Crop(basin, crop, df_summary):
    data = df_summary[(df_summary['Basin'] == basin) & (df_summary['Crop'] == crop)]
    if data.empty: return

    # 1. Fetch Boundaries
    n_boundary, p_boundary = Get_Boundaries(basin, crop)
    
    # 2. Logic for Scenarios
    scen3_n = data[data['Scenario'] == "Sus_Irri+Red_Fert"]['N_Runoff_ktons'].values
    show_fourth = (len(scen3_n) > 0 and n_boundary is not None and scen3_n[0] <= n_boundary)
    active_scenarios = SCENARIOS if show_fourth else SCENARIOS[:3]

    fig, axes = plt.subplots(1, 3, figsize=(8, 5)) 
    vars_to_plot = ['Production_ktons', 'N_Runoff_ktons', 'P_Runoff_ktons']
    
    # Corrected Boundaries list to match indices of vars_to_plot
    boundaries = [None, n_boundary, p_boundary]

    for i, var in enumerate(vars_to_plot):
        ax = axes[i]
        
        # 3. Unit Handling
        # Production -> Mt (div by 1000), Runoff -> ktons (div by 1)
        divisor = 1000 if var == 'Production_ktons' else 1
        
        vals = []
        for s in active_scenarios:
            v = data[data['Scenario'] == s][var].values
            vals.append(v[0]/divisor if len(v) > 0 else 0)
        
        base_val = vals[0]
        x_pos = np.arange(len(vals))
        
        # 4. Plot Bars
        for j in range(len(vals)):
            ax.bar(x_pos[j], vals[j], color=COLORS[j], 
                   edgecolor="white", linewidth=1.0, zorder=1, width=0.6)

        # 5. Percentage Labels (Conditional Decimals)
            if j > 0:
                diff = vals[j] - base_val
                pct = (diff / base_val * 100) if base_val != 0 else 0
                
                # Logic: If > 1% or < -1%, no decimal. If 0, "0.0%". Else, 1 decimal.
                if abs(pct) < 0.1:
                    label = "0%"
                elif abs(pct) >= 1.0:
                    label = f"{pct:+.0f}%" # No decimal
                else:
                    label = f"{pct:+.1f}%" # Keep decimal for small changes
                
                ax.text(j, vals[j] + (max(vals)*0.03), label, 
                        ha='center', va='bottom', fontsize=14, color='black')

        # 6. DRAW BOUNDARY LINE (Crucial Fix)
        # Ensure divisor is applied to the boundary value as well
        if boundaries[i] is not None:
            b_val = boundaries[i] / divisor
            # zorder=5 ensures it sits on top of everything
            ax.axhline(b_val, color="#b91515", linestyle='--', linewidth=2, zorder=5)

        # 7. Y-Axis Limits (Must include the boundary value!)
        plot_max = max(vals)
        if boundaries[i] is not None:
            plot_max = max(plot_max, boundaries[i]/divisor)
            
        upper_limit, tick_spacing = get_smart_ticks(plot_max)
        ax.set_ylim(0, upper_limit)
        ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))

        # 8. Aesthetic Cleanup
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks([]) 
        ax.tick_params(axis='y', labelsize=16)

    plt.tight_layout()
    plt.savefig(os.path.join(figure_dir, f"{basin}_{crop}_3SCE.png"), bbox_inches='tight')
    plt.close()


# --- Main Execution ---
# Load the summary CSV we created in the previous step
summary_file = os.path.join(Prod_Runoff_dir, "Production_Runoff_Summary.csv")
if os.path.exists(summary_file):
    df_summary = pd.read_csv(summary_file)

    for basin in STUDY_AREAS:
        for crop in all_crops:
            print(f"Plotting {basin} - {crop}...")
            Plot_Basin_Crop(basin, crop, df_summary)
    
    print(f"\nAll plots saved to: {figure_dir}")
else:
    print(f"Error: Summary file not found at {summary_file}")