# This code is used to plot
# 1. Total prouduction (sum of all crop types) (yaxis) and N (xaxis - left), P runoff (xaxis - right) changes when accepting different percentage of yield reduction 
# 2. Current boundaries and boundaries when excluding the contribution of sewage (dashed lines)

import os
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# --- 1. Configuration & Paths ---
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/5_TradeOffs"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Fig5"
current_boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/1_Boundary"
nosewage_boundary_dir = input_dir

os.makedirs(output_dir, exist_ok=True)

dissolved_N_frac = 0.7
dissolved_P_frac = 0.25

tradeoff_csv = os.path.join(input_dir, "Scenarios_Production_Runoff_Summary.csv")
Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]

CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}

SCENARIO_ORDER = {
    "Baseline": 0,
    # "10% Reduction": 1,
    "20% Reduction": 1,
    "30% Reduction": 2,
    "40% Reduction": 3,
    "50% Reduction": 4
}

# --- Helper Function for Dynamic Nice Numbers ---
def get_nice_axis_params(min_val, max_val):
    """
    Computes a clean step interval and nice bounds for any scale.
    """
    if max_val == min_val:
        return min_val, max_val + 1, 1
        
    span = max_val - min_val
    # Get the order of magnitude of the span
    exponent = math.floor(math.log10(span))
    fraction = span / (10 ** exponent)
    
    # Select a clean nice step based on fraction magnitude
    if fraction <= 1:
        step = 0.2
    elif fraction <= 2.5:
        step = 0.5
    elif fraction <= 5.0:
        step = 1.0
    else:
        step = 2.0
        
    step *= (10 ** exponent)
    
    # Round bounds to perfectly match the clean steps
    nice_min = math.floor(min_val / step) * step
    nice_max = math.ceil(max_val / step) * step
    
    return nice_min, nice_max, step


# --- 2. Load Scenario Data ---
if not os.path.exists(tradeoff_csv):
    raise FileNotFoundError(f"Missing scenario summary dataset at: {tradeoff_csv}")

df_all = pd.read_csv(tradeoff_csv)

# --- 3. Main Plotting Loop ---
for basin in Studyareas:
    print(f"Processing basin: {basin}...")
    
    df_basin = df_all[df_all['Basin'].str.lower() == basin.lower()].copy()
    if df_basin.empty:
        print(f"  [Warning] No scenario data found for {basin}. Skipping.")
        continue
        
    available_summary_crops = df_basin['Crop'].unique()
    df_basin_filtered = df_basin[df_basin['Crop'].isin(available_summary_crops)]
    
    df_agg = df_basin_filtered.groupby('Scenario')[['Production_ktons', 'N_Runoff_ktons', 'P_Runoff_ktons']].sum().reset_index()
    df_agg['sort_idx'] = df_agg['Scenario'].map(SCENARIO_ORDER)
    df_agg = df_agg.dropna(subset=['sort_idx']).sort_values('sort_idx')

    # Baseline scenario
    # baseline_series = df_agg.loc[df_agg['Scenario'] == 'Baseline', 'Production_ktons']
    # if not baseline_series.empty:
    #     baseline_prod = baseline_series.iloc[0]
    #     # Avoid division by zero if baseline production happens to be 0
    #     if baseline_prod != 0:
    #         prod_mt = 100 * df_agg['Production_ktons'] / baseline_prod
    #     else:
    #         prod_mt = df_agg['Production_ktons'] * 0.0
    # else:
    #     print(f"  [Warning] 'Baseline' scenario missing for basin {basin}. Falling back to default scaling.")
    #     prod_mt = 0.0
    
    # Unit Conversions: Only transform Crop Production (ktons -> Mt)
    prod_mt = df_agg['Production_ktons'] / 1000.0
    n_rf = df_agg['N_Runoff_ktons'] 
    p_rf = df_agg['P_Runoff_ktons'] 

    # --- 4. Boundary Loading ---
    basin_boundary_crops = []
    for s_crop in available_summary_crops:
        if s_crop in CROP_MAP:
            basin_boundary_crops.extend(CROP_MAP[s_crop])

    current_n_bound, current_p_bound = None, None
    nosewage_n_bound, nosewage_p_bound = None, None

    b_file_path = os.path.join(current_boundary_dir, f"{basin}_crop-specific_boundaries.csv")
    no_sewage_b_file_path = os.path.join(nosewage_boundary_dir, f"{basin}_crop-specific_boundaries.csv")   

    if os.path.exists(b_file_path) and os.path.exists(no_sewage_b_file_path):
        df_b = pd.read_csv(b_file_path, index_col=0)
        df_b_no_sewage = pd.read_csv(no_sewage_b_file_path, index_col=0)
        
        df_b.index = df_b.index.astype(str).str.strip()
        df_b_no_sewage.index = df_b_no_sewage.index.astype(str).str.strip()
        
        active_rows = [row for row in basin_boundary_crops if row in df_b.index]
        if active_rows:
            current_n_bound = dissolved_N_frac * float(df_b.loc[active_rows, 'N [ktons]'].sum())
            current_p_bound = dissolved_P_frac * float(df_b.loc[active_rows, 'P [ktons]'].sum())

        active_rows_ns = [row for row in basin_boundary_crops if row in df_b_no_sewage.index]
        if active_rows_ns:
            nosewage_n_bound = dissolved_N_frac * float(df_b_no_sewage.loc[active_rows_ns, 'N [ktons]'].sum())
            nosewage_p_bound = dissolved_P_frac * float(df_b_no_sewage.loc[active_rows_ns, 'P [ktons]'].sum())
    else:
        print(f"  [Error] Skipping {basin}: Boundary files missing from cluster drive paths.")
        continue

    # --- 5. Generate Stacked Layout Framework (N top, P bottom) ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10), sharex=True)
    color_n, color_p = '#CD2C58', '#F17727'
    
    # --- SUBPLOT 1 (Top): Nitrogen Runoff ---
    ax1.plot(prod_mt, n_rf, color=color_n, marker='o', markersize=15, linewidth=3)
    ax1.grid(True, linestyle=':', alpha=0.4)
    
    if current_n_bound and current_n_bound > 0:
        ax1.axhline(y=current_n_bound, color="black", linestyle='-', linewidth=2)
    if nosewage_n_bound and nosewage_n_bound > 0:
        ax1.axhline(y=nosewage_n_bound, color="black", linestyle='--', linewidth=2)
        
    # Calculate Nice Axis Bounds for Nitrogen
    n_all = list(n_rf) + [v for v in [current_n_bound, nosewage_n_bound] if v]
    n_min_nice, n_max_nice, n_step = get_nice_axis_params(min(n_all), max(n_all))
    
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
        
    # Calculate Nice Axis Bounds for Phosphorus
    p_all = list(p_rf) + [v for v in [current_p_bound, nosewage_p_bound] if v]
    p_min_nice, p_max_nice, p_step = get_nice_axis_params(min(p_all), max(p_all))
    
    ax2.set_ylim(p_min_nice, p_max_nice)
    ax2.yaxis.set_major_locator(MultipleLocator(p_step))
    
    # Calculate Nice Axis Bounds for Shared X-Axis (Production)
    x_min_nice, x_max_nice, x_step = get_nice_axis_params(min(prod_mt), max(prod_mt))
    # x_min_nice, x_max_nice, x_step = 50, 100, 10
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
    fig_name = f"TradeOff_Runoff_Boundaries_Stacked_{basin}.png"
    plt.savefig(os.path.join(output_dir, fig_name), dpi=300, bbox_inches='tight')
    plt.close()

print("All stacked minimalist charts exported successfully with clean step limits.")