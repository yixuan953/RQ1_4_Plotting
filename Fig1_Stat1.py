# This code is used to 
# 1. Summarize the boundaries for total and agricultural N and P load for the whole basin 
# 2. Output the results in .csv file


import os
import numpy as np
import xarray as xr 
import pandas as pd

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/2_Critical_NP_losses/Method3/Crit_Load"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/1_Boundary_load"
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def CalBasinSum(file_path):

    ds = xr.open_dataset(file_path)
    var_name = list(ds.data_vars)[0]
    
    # Select year 2005 and sum (handling potential dimension naming)
    data_2005 = ds[var_name].sel(year=2005)
    total_sum = data_2005.sum().values
    
    # Conversion: kg to ktons (10^6)
    value_ktons = total_sum * 1e-6
    
    ds.close()
    return value_ktons

for basin in Studyareas:
    print(f"Processing basin: {basin}...")
    
    # Define file paths
    total_N_file = os.path.join(input_dir, f"{basin}_total_crit_N_load.nc")
    total_P_file = os.path.join(input_dir, f"{basin}_total_crit_P_load.nc")
    agri_N_file = os.path.join(input_dir, f"{basin}_agri_crit_N_load.nc")
    agri_P_file = os.path.join(input_dir, f"{basin}_agri_crit_P_load.nc")

    # Calculate sums
    sum_total_N = CalBasinSum(total_N_file)
    sum_total_P = CalBasinSum(total_P_file)
    sum_agri_N = CalBasinSum(agri_N_file)
    sum_agri_P = CalBasinSum(agri_P_file)

    # Create a DataFrame for the CSV structure
    # Rows: Total, Agri | Columns: N, P
    data = {
        "N [ktons]": [sum_total_N, sum_agri_N],
        "P [ktons]": [sum_total_P, sum_agri_P]
    }
    df = pd.DataFrame(data, index=["Total", "Agri"])

    # Output to CSV
    output_file_name = os.path.join(output_dir, f"{basin}_crit_load_sum.csv")
    df.to_csv(output_file_name)
    
print("Processing complete. Files saved in:", output_dir)



