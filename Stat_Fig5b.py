# This script correctly sums up all crops within a basin FIRST, then calculates the basin-level metrics

import os
import numpy as np
import pandas as pd

# Define paths
input_file = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/3_TradeOffs/Scenarios_Production_Runoff_Summary.csv"
boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/1_Boundary"
output_file = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Reconcile/NPFootprint_End-of-Pipe_Removal_Ratios.csv"

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]

# Load the main simulated production and runoff data
df_sim = pd.read_csv(input_file)
df_sim["Basin_lower"] = df_sim["Basin"].str.lower()
df_sim["Crop_lower"] = df_sim["Crop"].str.lower()

# Map boundary file crop identifiers to simulation runoff crop names
crop_mapping = {
    "winterwheat": "Wheat",
    "maize": "Maize",
    "mainrice": "Rice",
    "secondrice": "Rice",
    "soybean": "Soybean",
}

# The explicit footprint reduction tiers requested for step 3
footprint_reductions = ["Current footprint", "Required footprint", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%"]
yield_scenarios = ["Baseline", "Sus_Irrigation", "5% Reduction", "10% Reduction", "20% Reduction", "30% Reduction", "40% Reduction", "50% Reduction"]

results = []

for studyarea in Studyarea:
    boundary_file = os.path.join(boundary_dir, f"{studyarea}_crop-specific_boundaries.csv")
    
    if not os.path.exists(boundary_file):
        print(f"Boundary file missing for {studyarea}, skipping...")
        continue
        
    df_bound = pd.read_csv(boundary_file)
    df_bound.columns = df_bound.columns.str.strip()
    id_col = df_bound.columns[0]
    df_bound[id_col] = df_bound[id_col].str.strip().str.lower()
    
    df_basin_sim = df_sim[df_sim["Basin_lower"] == studyarea.lower()]
    
    # --- STEP A: Sum up ALL crop boundaries for this basin first ---
    total_n_boundary = 0.0
    total_p_boundary = 0.0
    valid_sim_crops = []
    
    for bound_crop, sim_crop in crop_mapping.items():
        crop_bound_row = df_bound[df_bound[id_col] == bound_crop]
        if not crop_bound_row.empty:
            total_n_boundary += 0.7 * crop_bound_row["N [ktons]"].values[0]
            total_p_boundary += 0.25 * crop_bound_row["P [ktons]"].values[0]
            if sim_crop.lower() not in valid_sim_crops:
                valid_sim_crops.append(sim_crop.lower())
                
    # Filter basin data to only include these mapped crops
    df_basin_crops = df_basin_sim[df_basin_sim["Crop_lower"].isin(valid_sim_crops)].copy()
    
    # --- STEP B: Aggregate production & runoff across ALL crops for each scenario row ---
    df_aggregated = df_basin_crops.groupby("Scenario").agg({
        "Production_ktons": "sum",
        "N_Runoff_ktons": "sum",
        "P_Runoff_ktons": "sum"
    }).reset_index()
    
    # --- STEP C: Calculate metrics using the true Basin-wide totals ---
    for pct_yield in yield_scenarios:
        row = df_aggregated[df_aggregated["Scenario"] == pct_yield]
        
        if row.empty:
            continue
            
        prod = row["Production_ktons"].values[0]
        n_runoff = row["N_Runoff_ktons"].values[0]
        p_runoff = row["P_Runoff_ktons"].values[0]
        
        # 1) Base Footprint calculations for total basin
        n_footprint = (1000 * n_runoff / prod) if prod > 0 else 0.0
        p_footprint = (1000 * p_runoff / prod) if prod > 0 else 0.0
        
        # 2) Required Footprint metrics based on total basin boundaries
        n_req_footprint = (1000 * total_n_boundary / prod) if prod > 0 else 0.0
        p_req_footprint = (1000 * total_p_boundary / prod) if prod > 0 else 0.0
        
        # 3) Loop through each footprint reduction variant
        for ft_red in footprint_reductions:
            if ft_red == "Current footprint":
                n_decreased_ft = n_footprint
                p_decreased_ft = p_footprint
                scenario_label = f"{pct_yield} yield + Current footprint"
            elif ft_red == "Required footprint":
                n_decreased_ft = n_req_footprint
                p_decreased_ft = p_req_footprint
                scenario_label = f"{pct_yield} yield + Required footprint"
            else:
                factor = float(ft_red.replace("%", "")) / 100.0
                n_decreased_ft = n_footprint * (1.0 - factor)
                p_decreased_ft = p_footprint * (1.0 - factor)
                scenario_label = f"{pct_yield} yield + {ft_red} Footprint reduction"
            
            # End-of-pipe removal ratio calculation
            n_removal_ratio = ((n_decreased_ft - n_req_footprint) / n_decreased_ft) if n_decreased_ft > 0 else 0.0
            p_removal_ratio = ((p_decreased_ft - p_req_footprint) / p_decreased_ft) if p_decreased_ft > 0 else 0.0
            
            # Format to 0.00 if boundaries are already naturally satisfied
            n_removal_ratio = max(0.0, n_removal_ratio)
            p_removal_ratio = max(0.0, p_removal_ratio)
            
            results.append({
                "Basin": studyarea,
                "Scenario": scenario_label,
                "N footprint": round(n_decreased_ft, 3),
                "P footprint": round(p_decreased_ft, 3),
                "Required N footprint": round(n_req_footprint, 3),
                "Required P footprint": round(p_req_footprint, 3),
                "N Boundary removal ratio": round(n_removal_ratio, 4),
                "P Boundary removal ratio": round(p_removal_ratio, 4)
            })

# Save final structured outputs
if len(results) > 0:
    output_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    output_df.to_csv(output_file, index=False)
    print(f"Data mapping completed. Output saved to: {output_file}")