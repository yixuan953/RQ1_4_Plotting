# This code is used to 
# 1 - Summarize what is the boundaries for irrigation, N and P runoff for each crop types (dissolvable + particle)
# 2 - Input the results in the .csv from the previous step

import pandas as pd
import numpy as np
import xarray as xr
import os

# Paths for data input (NetCDFs)
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

# Paths for CSV processing
csv_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/4_Reconcile/75Sewage"

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

prop_dissolved_N = 0.7
prop_dissolved_P = 0.25

def GetCropWNP(file_name, low_runoff_path):
    with xr.open_dataset(file_name) as ds:
        # Use .values to ensure we get a clean sum without coord issues
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
        
        with xr.open_dataset(low_runoff_path) as ds_lr:
            low_runoff = ds_lr["Low_Runoff"]
            mask = mask.where(low_runoff != 1, np.nan)
        
        # Calculate load and convert kg to ktons (1e-6)
        CWATER = (ds["Sus_irrigation_amount"] * mask).sum().values
        CN = (ds["Crit_N_Runoff"] * mask).sum().values * 1e-6 * (1/prop_dissolved_N)
        CP = (ds["Crit_P_Runoff"] * mask).sum().values * 1e-6 * (1/prop_dissolved_P)

    return float(CN), float(CP), float(CWATER)

for basin in Studyarea:
    print(f"Processing basin: {basin}")
    
    # 1. Path to your existing CSV (the one with Total and Agri)
    existing_csv_path = os.path.join(csv_dir, f"{basin}_boundaries_sum.csv")
    
    if not os.path.exists(existing_csv_path):
        print(f"Warning: Existing CSV not found for {basin}, skipping...")
        continue

    # 2. Load existing CSV
    # We assume first column is the index (Total, Agri)
    df_existing = pd.read_csv(existing_csv_path, index_col=0)
    
    # 3. Calculate new crop rows
    crop_rows = []
    for crop in Croptypes:
        file_name = os.path.join(input_dir, f"{basin}_{crop}_summary.nc")
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
        
        if os.path.exists(file_name) and os.path.exists(low_runoff_path):
            val_n, val_p, val_water = GetCropWNP(file_name, low_runoff_path)
            
            crop_rows.append({
                "Category": crop, # This matches the index style
                "Water [m³]": val_water,
                "N [ktons]": val_n,
                "P [ktons]": val_p
            })
    
    # 4. Convert new crops to DataFrame and Append
    if crop_rows:
        df_new_crops = pd.DataFrame(crop_rows).set_index("Category")
        
        # FIX: Print columns to debug if a mismatch happens again
        if len(df_new_crops.columns) != len(df_existing.columns):
            print(f"⚠️ Column count mismatch for {basin}!")
            print(f"Existing columns ({len(df_existing.columns)}): {list(df_existing.columns)}")
            print(f"New columns ({len(df_new_crops.columns)}): {list(df_new_crops.columns)}")
            print("Attempting to align columns based on existing names...")
            
            # Robust fix: Map your new data directly to whatever the existing column names are
            # This assumes your order is Water, N, P in both files
            if len(df_existing.columns) == 3:
                df_new_crops.columns = df_existing.columns
            else:
                # If the existing CSV genuinely only has 2 data columns, 
                # we match by name rather than forcing them.
                pass 
        else:
            # If lengths match perfectly, safely assign
            df_new_crops.columns = df_existing.columns
        
        # Concatenate safely (axis=0 means stacking rows)
        # using sort=False keeps the original column order of df_existing
        df_final = pd.concat([df_existing, df_new_crops], axis=0, sort=False)
        
        # 5. Save back to CSV
        output_path = os.path.join(csv_dir, f"{basin}_crop-specific_boundaries.csv")
        df_final.to_csv(output_path)
        print(f"Updated CSV saved to: {output_path}")

print("Done!")


