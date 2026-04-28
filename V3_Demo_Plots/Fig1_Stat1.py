# This code is used to 
# 1. Summarize the boundaries for total and agricultural N and P load for the whole basin 
# 2. Output the results in .csv file


import os
import numpy as np
import xarray as xr 
import pandas as pd

model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/2_Critical_NP_losses/Method3/Crit_Load"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/1_Boundary_load"
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def SelYearCalBasinSum(file_path, mask):

    ds = xr.open_dataset(file_path)
    var_name = list(ds.data_vars)[0]
    
    # Select year 2015 and sum (handling potential dimension naming)
    data_2015 = ds[var_name].sel(year=2014)
    data_2015_mask = mask * data_2015
    total_sum = data_2015_mask.sum().values
    
    # Conversion: kg to ktons (10^6)
    value_ktons = total_sum * 1e-6
    
    ds.close()
    return value_ktons

def CalBasinSum(file_path, mask):

    ds = xr.open_dataset(file_path)
    var_name = list(ds.data_vars)[0]
    
    # Select year 2015 and sum (handling potential dimension naming)
    data_2015 = ds[var_name]
    data_2015_mask = mask * data_2015
    total_sum = data_2015_mask.sum().values
    
    # Conversion: kg to ktons (10^6)
    value_ktons = total_sum * 1e-6
    
    ds.close()
    return value_ktons

for basin in Studyareas:
    print(f"Processing basin: {basin}...")

    file_name = os.path.join(model_summary_dir, f"{basin}_maize_summary.nc")
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")

    with xr.open_dataset(file_name) as ds:
        # Use .values to ensure we get a clean sum without coord issues
        basin_mask = ds["Basin_mask"]
        
    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]
        mask = basin_mask.where(low_runoff != 1, np.nan)

    # Define file paths
    total_N_file = os.path.join(input_dir, f"{basin}_total_crit_N_load.nc")
    total_P_file = os.path.join(input_dir, f"{basin}_total_crit_P_load.nc")
    agri_N_file = os.path.join(input_dir, f"{basin}_agri_crit_N_load.nc")
    agri_P_file = os.path.join(input_dir, f"{basin}_agri_crit_P_load.nc")
    crop_N_file = os.path.join(input_dir, f"{basin}_cropland_crit_N_load.nc")
    crop_P_file = os.path.join(input_dir, f"{basin}_cropland_crit_P_load.nc")

    # Calculate sums
    sum_total_N = SelYearCalBasinSum(total_N_file, basin_mask)
    sum_total_P = SelYearCalBasinSum(total_P_file, basin_mask)
    sum_agri_N = SelYearCalBasinSum(agri_N_file, basin_mask)
    sum_agri_P = SelYearCalBasinSum(agri_P_file, basin_mask)
    sum_crop_N = CalBasinSum(crop_N_file, mask)
    sum_crop_P = CalBasinSum(crop_P_file, mask)

    # Create a DataFrame for the CSV structure
    # Rows: Total, Agri | Columns: N, P
    data = {
        "N [ktons]": [sum_total_N, sum_agri_N, sum_crop_N],
        "P [ktons]": [sum_total_P, sum_agri_P, sum_crop_P*4]
    }
    df = pd.DataFrame(data, index=["All sources", "Agriculture", "Cropland"])

    # Output to CSV
    output_file_name = os.path.join(output_dir, f"{basin}_crit_load_sum.csv")
    df.to_csv(output_file_name)
    
print("Processing complete. Files saved in:", output_dir)
