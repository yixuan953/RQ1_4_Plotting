# This code is used to plot
# 1. Total prouduction and N, P runoff changes when accepting different percentage of yield reduction 
# 2. Current boundaries and boundaries when excluding the contribution of sewage (dashed lines)

import os
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# --- 1. Configuration & Paths ---
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/3_TradeOffs"
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
    "Sus_Irrigation": 1,
    "5% Reduction": 2,
    "10% Reduction": 3,
    "20% Reduction": 4,
    "30% Reduction": 5,
    "40% Reduction": 6,
    "50% Reduction": 7
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
    
    df_agg_raw = df_basin_filtered.groupby('Scenario')[['Production_ktons', 'N_Runoff_ktons', 'P_Runoff_ktons']].sum().reset_index()
  
    Yp_series = df_agg_raw.loc[df_agg_raw['Scenario'] == 'Yp', 'Production_ktons']
    if not Yp_series.empty:
        baseline_prod = Yp_series.iloc[0]
    else:
        baseline_series = df_agg_raw.loc[df_agg_raw['Scenario'] == 'Baseline', 'Production_ktons']
        if not baseline_series.empty:
            baseline_prod = baseline_series.iloc[0]
            print(f"  [Warning] 'Yp' scenario missing for basin {basin}. Using 'Baseline' as fallback.")
        else:
            baseline_prod = 1.0
            print(f"  [Warning] Neither 'Yp' nor 'Baseline' found for basin {basin}. Scaling factor set to 1.0.")
    
    df_agg = df_agg_raw.copy()
    df_agg['sort_idx'] = df_agg['Scenario'].map(SCENARIO_ORDER)
    df_agg = df_agg.dropna(subset=['sort_idx']).sort_values('sort_idx')

    if baseline_prod != 0:
        prod_pct = 100.0 * df_agg['Production_ktons'] / baseline_prod
    else:
        prod_pct = df_agg['Production_ktons'] * 0.0
    
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
    
    # Dynamic Buffer calculation for centering target elements nicely
    y_buffer_ratio = 0.25 
    x_buffer_ratio = 0.08  # Increased slightly to prevent large marker clipping
    
    # Define custom colors for the first two specific scenario points
    pt_colors_n = ['grey' if i == 0 else '#0E75A9' if i == 1 else color_n for i in range(len(prod_pct))]
    pt_colors_p = ['grey' if i == 0 else '#0E75A9' if i == 1 else color_p for i in range(len(prod_pct))]

    # --- SUBPLOT 1 (Top): Nitrogen Runoff ---
    ax1.plot(prod_pct, n_rf, color="grey", linestyle='--',linewidth=3, zorder=4)
    ax1.scatter(prod_pct, n_rf, color=pt_colors_n, s=240, zorder=5)
    ax1.grid(True, linestyle=':', alpha=0.3, zorder=1)
    
    # Calculate Dynamic Limits for Y-axis
    n_all = list(n_rf) + [v for v in [current_n_bound, nosewage_n_bound] if v]
    n_min_raw, n_max_raw = min(n_all), max(n_all)
    n_span = n_max_raw - n_min_raw if n_max_raw != n_min_raw else 1.0
    
    n_min_nice, n_max_nice, n_step = get_nice_axis_params(
        n_min_raw - (n_span * y_buffer_ratio), 
        n_max_raw + (n_span * y_buffer_ratio)
    )
    ax1.set_ylim(n_min_nice, n_max_nice)
    
    # Shadows & Threshold Lines (Using a lighter, cleaner alpha and subtle hatch)
    if nosewage_n_bound and nosewage_n_bound > 0:
        ax1.axhspan(ymin=nosewage_n_bound, ymax=n_max_nice, color='#FF4D4D', alpha=0.08, hatch='//', zorder=2)
        ax1.axhline(y=nosewage_n_bound, color="red", linestyle='--', linewidth=2, zorder=3)
        
    if current_n_bound and current_n_bound > 0:
        ax1.axhline(y=current_n_bound, color="red", linestyle='-', linewidth=2, zorder=3)
        
    ax1.yaxis.set_major_locator(MultipleLocator(n_step))
    ax1.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
    ax1.tick_params(axis='y', labelsize=25)
    
    # --- SUBPLOT 2 (Bottom): Phosphorus Runoff ---
    ax2.plot(prod_pct, p_rf, color="grey", linestyle='--', linewidth=3, zorder=4)
    ax2.scatter(prod_pct, p_rf, color=pt_colors_p, s=240, zorder=5) 
    ax2.grid(True, linestyle=':', alpha=0.3, zorder=1)
    
    # Calculate Dynamic Limits for Y-axis
    p_all = list(p_rf) + [v for v in [current_p_bound, nosewage_p_bound] if v]
    p_min_raw, p_max_raw = min(p_all), max(p_all)
    p_span = p_max_raw - p_min_raw if p_max_raw != p_min_raw else 1.0
    
    p_min_nice, p_max_nice, p_step = get_nice_axis_params(
        p_min_raw - (p_span * y_buffer_ratio), 
        p_max_raw + (p_span * y_buffer_ratio)
    )
    ax2.set_ylim(p_min_nice, p_max_nice)
    
    # Shadows & Threshold Lines
    if nosewage_p_bound and nosewage_p_bound > 0:
        ax2.axhspan(ymin=nosewage_p_bound, ymax=p_max_nice, color='#FF4D4D', alpha=0.08, hatch='//', zorder=2)
        ax2.axhline(y=nosewage_p_bound, color="red", linestyle='--', linewidth=2, zorder=3)
        
    if current_p_bound and current_p_bound > 0:
        ax2.axhline(y=current_p_bound, color="red", linestyle='-', linewidth=2, zorder=3)
        
    ax2.yaxis.set_major_locator(MultipleLocator(p_step))
    
    # --- Global X-axis adjustments ---
    x_min, x_max = min(prod_pct), max(prod_pct)
    x_span = x_max - x_min if x_max != x_min else 1.0
    x_min_nice, x_max_nice, x_step = get_nice_axis_params(
        x_min - (x_span * x_buffer_ratio),
        x_max + (x_span * x_buffer_ratio)
    )
    ax2.set_xlim(x_min_nice, x_max_nice)
    ax2.xaxis.set_major_locator(MultipleLocator(x_step))
    
    # Invert x-axis so Baseline (high production) is on the right
    ax2.invert_xaxis()
    
    # ADDED 'pad' HERE: This pushes the numbers slightly away from the axis line 
    # to stop the corner values from colliding into each other.
    ax2.tick_params(axis='x', which='both', bottom=True, labelbottom=True, labelsize=25, pad=10)
    ax1.tick_params(axis='y', labelsize=25, pad=8)
    ax2.tick_params(axis='y', labelsize=25, pad=8)

    # --- 6. Clean Frame Cleanup ---
    for ax in [ax1, ax2]:
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel("")
        for spine in ax.spines.values():
            spine.set_color('#333333')
            spine.set_linewidth(1.2)

    # Automatically realign the subplot labels dynamically before exporting
    fig.align_labels()

    # --- 7. Export ---
    # Change tight_layout padding slightly if needed, or stick with your setup:
    plt.tight_layout(pad=1.5)
    fig_name = f"TradeOff_Runoff_Boundaries_Stacked_{basin}.png"
    plt.savefig(os.path.join(output_dir, fig_name), dpi=300, bbox_inches='tight')
    plt.close()

print("All stacked minimalist charts exported successfully with clean step limits.")