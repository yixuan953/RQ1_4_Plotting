# This code is used to plot
# 1. Production (yaxis) and N (xaxis - left), P runoff (xaxis - right) changes per crop type
# 2. Current boundaries and boundaries when excluding the contribution of sewage (dashed lines)

import os
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# --- 1. Configuration & Paths ---
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/5_TradeOffs"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Extended_Fig/Fig6"
current_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/1_Boundary"
nosewage_boundary_dir = input_dir

os.makedirs(output_dir, exist_ok=True)

dissolved_N_frac = 0.7
dissolved_P_frac = 0.25

tradeoff_csv = os.path.join(input_dir, "Scenarios_Production_Runoff_Summary.csv")
Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]

TARGET_CROPS = ["Wheat", "Maize", "Rice", "Soybean"]

CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}

SCENARIO_ORDER = {
    "Baseline": 0,
    "20% Reduction": 1,
    "30% Reduction": 2,
    "40% Reduction": 3,
    "50% Reduction": 4
}

# --- Helper Function for Clean Axis with Max 4 Ticks ---
def get_nice_axis_params(min_val, max_val, max_ticks=4):
    """
    Computes a clean step interval and nice bounds, strictly limiting
    the number of major tick intervals to `max_ticks` (default 4).
    """
    if max_val == min_val:
        return min_val, max_val + 1, 1
        
    span = max_val - min_val
    exponent = math.floor(math.log10(span))
    fraction = span / (10 ** exponent)
    
    # Establish a starting base step based on standard graph intervals
    if fraction <= 1.5:
        step_base = 0.5
    elif fraction <= 3.0:
        step_base = 1.0
    elif fraction <= 6.0:
        step_base = 2.0
    else:
        step_base = 5.0
        
    step = step_base * (10 ** exponent)
    
    # Calculate perfect rounded bounds based on initial step choice
    nice_min = math.floor(min_val / step) * step
    nice_max = math.ceil(max_val / step) * step
    
    # Check if this results in too many ticks. If so, upscale the interval.
    while ((nice_max - nice_min) / step) >= max_ticks:
        if step_base == 0.5:
            step_base = 1.0
        elif step_base == 1.0:
            step_base = 2.0
        elif step_base == 2.0:
            step_base = 5.0
        else:
            step_base *= 2.0  # Keep expanding by factors of 2 if data is extremely skewed
            
        step = step_base * (10 ** exponent)
        nice_min = math.floor(min_val / step) * step
        nice_max = math.ceil(max_val / step) * step

    return nice_min, nice_max, step


# --- 2. Load Scenario Data ---
if not os.path.exists(tradeoff_csv):
    raise FileNotFoundError(f"Missing scenario summary dataset at: {tradeoff_csv}")

df_all = pd.read_csv(tradeoff_csv)

# --- 3. Main Plotting Loops ---
for basin in Studyareas:
    print(f"Processing basin: {basin}...")
    
    df_basin = df_all[df_all['Basin'].str.lower() == basin.lower()].copy()
    if df_basin.empty:
        print(f"  [Warning] No scenario data found for {basin}. Skipping.")
        continue
        
    for crop_name in TARGET_CROPS:
        df_crop = df_basin[df_basin['Crop'] == crop_name].copy()
        
        if df_crop.empty:
            continue
            
        df_agg = df_crop.groupby('Scenario')[['Production_ktons', 'N_Runoff_ktons', 'P_Runoff_ktons']].sum().reset_index()
        df_agg['sort_idx'] = df_agg['Scenario'].map(SCENARIO_ORDER)
        df_agg = df_agg.dropna(subset=['sort_idx']).sort_values('sort_idx')
        
        if df_agg.empty:
            continue

        prod_mt = df_agg['Production_ktons'] / 1000.0
        n_rf = df_agg['N_Runoff_ktons'] 
        p_rf = df_agg['P_Runoff_ktons'] 

        # --- 4. Crop-Specific Boundary Loading ---
        current_n_bound, current_p_bound = None, None
        nosewage_n_bound, nosewage_p_bound = None, None
        
        sub_crops = CROP_MAP[crop_name]

        b_file_path = os.path.join(current_boundary_dir, f"{basin}_crop-specific_boundaries.csv")
        no_sewage_b_file_path = os.path.join(nosewage_boundary_dir, f"{basin}_crop-specific_boundaries.csv")   

        if os.path.exists(b_file_path) and os.path.exists(no_sewage_b_file_path):
            df_b = pd.read_csv(b_file_path, index_col=0)
            df_b_no_sewage = pd.read_csv(no_sewage_b_file_path, index_col=0)
            
            df_b.index = df_b.index.astype(str).str.strip()
            df_b_no_sewage.index = df_b_no_sewage.index.astype(str).str.strip()
            
            active_rows = [row for row in sub_crops if row in df_b.index]
            if active_rows:
                current_n_bound = dissolved_N_frac * float(df_b.loc[active_rows, 'N [ktons]'].sum())
                current_p_bound = dissolved_P_frac * float(df_b.loc[active_rows, 'P [ktons]'].sum())

            active_rows_ns = [row for row in sub_crops if row in df_b_no_sewage.index]
            if active_rows_ns:
                nosewage_n_bound = dissolved_N_frac * float(df_b_no_sewage.loc[active_rows_ns, 'N [ktons]'].sum())
                nosewage_p_bound = dissolved_P_frac * float(df_b_no_sewage.loc[active_rows_ns, 'P [ktons]'].sum())

        # --- 5. Generate Stacked Layout Framework ---
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10), sharex=True)
        color_n, color_p = '#CD2C58', '#F17727'
        
        # --- SUBPLOT 1 (Top): Nitrogen Runoff ---
        ax1.plot(prod_mt, n_rf, color=color_n, marker='o', markersize=15, linewidth=3)
        ax1.grid(True, linestyle=':', alpha=0.4)
        
        if current_n_bound and current_n_bound > 0:
            ax1.axhline(y=current_n_bound, color="black", linestyle='-', linewidth=2)
        if nosewage_n_bound and nosewage_n_bound > 0:
            ax1.axhline(y=nosewage_n_bound, color="black", linestyle='--', linewidth=2)
            
        # Dynamically scale N axis down to a max of 4 ticks
        n_all = list(n_rf) + [v for v in [current_n_bound, nosewage_n_bound] if v and v > 0]
        n_min_nice, n_max_nice, n_step = get_nice_axis_params(min(n_all), max(n_all), max_ticks=4)
        
        ax1.set_ylim(n_min_nice, n_max_nice)
        ax1.yaxis.set_major_locator(MultipleLocator(n_step))
        ax1.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
        ax1.tick_params(axis='y', labelsize=25)
        

        # --- SUBPLOT 2 (Bottom): Phosphorus Runoff ---
        ax2.plot(prod_mt, p_rf, color=color_p, marker='o', markersize=15, linewidth=3)
        ax2.grid(True, linestyle=':', alpha=0.4)
        
        if current_p_bound and current_p_bound > 0:
            ax2.axhline(y=current_p_bound, color="black", linestyle='-', linewidth=2)
        if nosewage_p_bound and nosewage_p_bound > 0:
            ax2.axhline(y=nosewage_p_bound, color="black", linestyle='--', linewidth=2)
            
        # Dynamically scale P axis down to a max of 4 ticks
        p_all = list(p_rf) + [v for v in [current_p_bound, nosewage_p_bound] if v and v > 0]
        p_min_nice, p_max_nice, p_step = get_nice_axis_params(min(p_all), max(p_all), max_ticks=4)
        
        ax2.set_ylim(p_min_nice, p_max_nice)
        ax2.yaxis.set_major_locator(MultipleLocator(p_step))
        
        # Shared X-Axis (Production) down to a max of 4 ticks
        x_min_nice, x_max_nice, x_step = get_nice_axis_params(min(prod_mt), max(prod_mt), max_ticks=4)
        ax2.set_xlim(x_min_nice, x_max_nice)
        ax2.xaxis.set_major_locator(MultipleLocator(x_step))
        
        ax2.tick_params(axis='x', which='both', bottom=True, labelbottom=True, labelsize=25)
        ax2.tick_params(axis='y', labelsize=25)

        # --- 6. Complete Removal of Titles, Text labels, and Legends ---
        ax1.set_title("")
        ax1.set_xlabel("")
        ax1.set_ylabel("")
        ax2.set_title("")
        ax2.set_xlabel("")
        ax2.set_ylabel("")

        # --- 7. Tight Clean Export ---
        plt.tight_layout()
        fig_name = f"TradeOff_Runoff_Boundaries_Stacked_{basin}_{crop_name}.png"
        plt.savefig(os.path.join(output_dir, fig_name), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    Saved cleanly scaled chart for {crop_name} in {basin}")

print("All crop-specific stacked charts exported successfully.")