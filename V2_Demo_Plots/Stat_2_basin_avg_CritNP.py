# This script is used to calculate the average critical N, P runoff [kg/ha] for each basin and each croptype, and save it in to .csv file

import pandas as pd
import numpy as np
import xarray as xr
import os

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V1_Statistics/2_Crit_NP_runoff"

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
    CN = ds["Crit_N_Runoff"] * mask
    CP = ds["Crit_P_Runoff"] * mask

    Basin_avg_CN = float(CN.sum()) / float(Total_HA.sum())
    Basin_avg_CP = float(CP.sum()) / float(Total_HA.sum())

    return Basin_avg_CN, Basin_avg_CP

for basin in Studyarea:
    
    basin_results = []
    
    for crop in Croptypes:
        file_name = os.path.join(input_dir, f"{basin}_{crop}_summary.nc")
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
        
        if os.path.exists(file_name) and os.path.exists(low_runoff_path):
            Basin_avg_CN, Basin_avg_CP = GetCropWNP(file_name, low_runoff_path)
            
            row = {
                "Crop type": crop,
                "Basin avg Crit N Runoff (kg/ha)": Basin_avg_CN,
                "Basin avg Crit P Runoff (kg/ha)": Basin_avg_CP
            }
            basin_results.append(row)
    
    if basin_results:
        df_basin = pd.DataFrame(basin_results)
        output_file = os.path.join(output_dir, f"{basin}_avg_Crit_NP_Runoff.csv")
        df_basin.to_csv(output_file, index=False)
        print(f"Saved average critical N and P runoff for {basin} to {output_file}")