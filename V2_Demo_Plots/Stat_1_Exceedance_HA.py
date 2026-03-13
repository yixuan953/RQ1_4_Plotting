# This script is used to calculate how many ha of harvested area has exceeded each boundary (irrigation, N runoff, P runoff) for each crop and each basin

import pandas as pd
import numpy as np
import xarray as xr
import os

model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V1_Statistics/3_Exceedance_HA_after_reduction"

os.makedirs(output_dir, exist_ok=True)

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean", "all_crop"]

def GetCropWNP(file_name, low_runoff_path):

    with xr.open_dataset(file_name) as ds:
        mask = xr.where(ds["Total_HA"] > 2500, 1, np.nan)

    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]
        
    mask = mask.where(low_runoff != 1, np.nan)
        
    # Extracting variables
    Total_HA = ds["Total_HA"] * mask
    Irri_HA = ds["Irrigated_HA"] * mask
    T_Irr = ds["Total_irrigation_amount"] * mask
    S_Irr = ds["Sus_irrigation_amount"] * mask
    N = ds["N_Runoff"] * mask
    CN = ds["Crit_N_Runoff"] * mask
    P = ds["P_Runoff"] * mask
    CP = ds["Crit_P_Runoff"] * mask

    return Total_HA, Irri_HA, T_Irr, S_Irr, N, P, CN, CP

# Loop through each basin to create individual CSVs
for basin in Studyarea:

    basin_results = []
    
    for crop in Croptypes:
        file_name = os.path.join(model_summary_dir, f"{basin}_{crop}_summary.nc")
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
        
        if os.path.exists(file_name) and os.path.exists(low_runoff_path):
            HA, Irri_HA, T_Irr, S_Irr, N, P, CN, CP = GetCropWNP(file_name, low_runoff_path)
            
            # Boolean masks
            w_exc = (T_Irr > S_Irr)
            n_exc = (N > CN)
            p_exc = (P > CP)

            total_ha_sum = float(HA.sum())
            
            if total_ha_sum > 0:
                row = {
                    "Crop type": crop,
                    "Irri": float(Irri_HA.where(w_exc).sum()) / total_ha_sum,
                    "N": float(HA.where(n_exc).sum()) / total_ha_sum,
                    "P": float(HA.where(p_exc).sum()) / total_ha_sum,
                    "Irri & N": float(Irri_HA.where(w_exc & n_exc).sum()) / total_ha_sum,
                    "Irri & P": float(Irri_HA.where(w_exc & p_exc).sum()) / total_ha_sum,
                    "N & P": float(HA.where(n_exc & p_exc).sum()) / total_ha_sum,
                    "All three": float(Irri_HA.where(w_exc & n_exc & p_exc).sum()) / total_ha_sum
                }
                basin_results.append(row)

    # Process the results for this specific basin
    if basin_results:
        df_basin = pd.DataFrame(basin_results)
        
        # Save unique file for this basin
        output_file = os.path.join(output_dir, f"{basin}_boundaries_ExceedanceHA.csv")
        df_basin.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")