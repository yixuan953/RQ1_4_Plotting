# This code is used to 
# 1 - Summarize what is the boundaries for N and P runoff for each crop types (dissolvable)
# 2 - Input the results in the .csv from the previous step

import pandas as pd
import numpy as np
import xarray as xr
import os

# Paths for data input (NetCDFs)
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

# Paths for CSV processing
csv_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/1_Boundary_load"

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

def GetCropWNP(file_name, low_runoff_path):
    with xr.open_dataset(file_name) as ds:
        # Use .values to ensure we get a clean sum without coord issues
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
        
        with xr.open_dataset(low_runoff_path) as ds_lr:
            low_runoff = ds_lr["Low_Runoff"]
            mask = mask.where(low_runoff != 1, np.nan)
        
        # Calculate load and convert kg to ktons (1e-6)
        CN = (ds["Crit_N_Runoff"] * mask).sum().values * 1e-6
        CP = (ds["Crit_P_Runoff"] * mask * 4).sum().values * 1e-6

    return float(CN), float(CP)

for basin in Studyarea:
    print(f"Processing basin: {basin}")
    
    # 1. Path to your existing CSV (the one with Total and Agri)
    existing_csv_path = os.path.join(csv_dir, f"{basin}_crit_load_sum.csv")
    
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
            val_n, val_p = GetCropWNP(file_name, low_runoff_path)
            
            crop_rows.append({
                "Category": crop, # This matches the index style
                "N [ktons]": val_n,
                "P [ktons]": val_p
            })
    
    # 4. Convert new crops to DataFrame and Append
    if crop_rows:
        df_new_crops = pd.DataFrame(crop_rows).set_index("Category")
        
        # Ensure column names match exactly for a clean append
        df_new_crops.columns = df_existing.columns
        
        # Concatenate existing and new
        df_final = pd.concat([df_existing, df_new_crops])
        
        # 5. Save back to CSV (overwriting or saving as a new version)
        output_path = os.path.join(csv_dir, f"{basin}_crit_load_sum_updated.csv")
        df_final.to_csv(output_path)
        print(f"Updated CSV saved to: {output_path}")

print("Done!")


