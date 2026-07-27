# This script is used to calculate the N, P delivery from difference sources in the two basins

import pandas as pd
import numpy as np
import xarray as xr
import os

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
# Outpu_directory
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/SM/2_NP_Delivery_IMAGE"

# Input directory for N
IMAGE_file_N = "/lustre/nobackup/WUR/ESG/zhou111/Data/Raw/IMAGE/Output-IMAGE_GNM-SSP1_oct2020-Nitrogen_Rivers-v2.nc"
ds_IMAGE_N = xr.open_dataset(IMAGE_file_N)
Agri_Surface_N = ds_IMAGE_N["Nsurface_runoff_agri"] 
Agri_Groundwater_N = ds_IMAGE_N["Ngroundwater_agri"]  
Sewage_N = ds_IMAGE_N["Nsewage"] 
Aquaculture_N = ds_IMAGE_N["Naquaculture"] 
total_N_load_non_manageable = ds_IMAGE_N["Ndeposition_water"]  + ds_IMAGE_N["Ngroundwater_nat"] + ds_IMAGE_N["Nsurface_runoff_nat"] + ds_IMAGE_N["Nvegetation"] 

total_N_load_non_manageable = total_N_load_non_manageable.sel(time = "2015-05-01")
Agri_Surface_N = Agri_Surface_N.sel(time = "2015-05-01")
Agri_Groundwater_N = Agri_Groundwater_N.sel(time = "2015-05-01")
Sewage_N = Sewage_N.sel(time = "2015-05-01")
Aquaculture_N = Aquaculture_N.sel(time = "2015-05-01")

lat = ds_IMAGE_N["lat"].values
lon = ds_IMAGE_N["lon"].values

# Input directory for P
IMAGE_file_P = "/lustre/nobackup/WUR/ESG/zhou111/Data/Raw/IMAGE/Output-IMAGE_GNM-SSP1_oct2020-Phosphorus_Rivers-v2.nc"
ds_IMAGE_P = xr.open_dataset(IMAGE_file_P)
total_P_load_non_manageable = ds_IMAGE_P["Psurface_runoff_nat"]  + ds_IMAGE_P["Pvegetation"]  + ds_IMAGE_P["Pweathering"]
Agri_Surface_P = ds_IMAGE_P["Psurface_runoff_agri"]
Sewage_P = ds_IMAGE_P["Psewage"]
Aquaculture_P = ds_IMAGE_P["Paquaculture"]

total_P_load_non_manageable = total_P_load_non_manageable.sel(time = "2015-05-01")
Agri_Surface_P = Agri_Surface_P.sel(time = "2015-05-01")
Sewage_P = Sewage_P.sel(time = "2015-05-01")
Aquaculture_P = Aquaculture_P.sel(time = "2015-05-01")

def CalBasinSum(data_2015, mask):

    data_2015_mask = mask * data_2015
    total_sum = data_2015_mask.sum().values
    
    # Conversion: kg to ktons (10^6)
    value_ktons = total_sum * 1e-6
    
    ds.close()
    return value_ktons

for basin in Studyarea:
    print(f"Processing basin: {basin}...")

    file_name = os.path.join(model_summary_dir, f"{basin}_maize_summary.nc")

    with xr.open_dataset(file_name) as ds:
        # Use .values to ensure we get a clean sum without coord issues
        mask = ds["Basin_mask"]

    sum_Agri_Surface_N = CalBasinSum(Agri_Surface_N, mask)
    sum_Agri_Groundwater_N = CalBasinSum(Agri_Groundwater_N, mask)
    sum_Sewage_N = CalBasinSum(Sewage_N, mask)
    sum_Aquaculture_N = CalBasinSum(Aquaculture_N, mask)
    sum_total_N_non_manageable = CalBasinSum(total_N_load_non_manageable, mask)

    sum_Agri_Surface_P = CalBasinSum(Agri_Surface_P, mask)
    sum_Agri_Groundwater_P = 0  # Assuming no groundwater P in the dataset
    sum_Sewage_P = CalBasinSum(Sewage_P, mask)
    sum_Aquaculture_P = CalBasinSum(Aquaculture_P, mask)
    sum_total_P_non_manageable = CalBasinSum(total_P_load_non_manageable, mask)

    # Create a DataFrame for the CSV structure
    # Rows: Total, Agri | Columns: N, P
    data = {
        "N [ktons]": [sum_total_N_non_manageable, sum_Agri_Surface_N, sum_Agri_Groundwater_N, sum_Sewage_N, sum_Aquaculture_N],
        "P [ktons]": [sum_total_P_non_manageable, sum_Agri_Surface_P, sum_Agri_Groundwater_P, sum_Sewage_P, sum_Aquaculture_P]
    }
    df = pd.DataFrame(data, index=["Non manageable sources", "Agri_Runoff", "Agri_Groundwater", "Sewage", "Aquaculture"])

    # Output to CSV
    output_file_name = os.path.join(output_dir, f"{basin}_IMAGE_NP_delivery_sum.csv")
    df.to_csv(output_file_name)
    
print("Processing complete. Files saved in:", output_dir)
