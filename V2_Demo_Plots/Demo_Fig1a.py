# This script is used to output the area where Water, N, P boundaries are exceeded (4 major crops) and export the .nc file
import os
import xarray as xr
import numpy as np

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

model_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
output_dir = model_output_dir

os.makedirs(output_dir, exist_ok=True)

def GetCropWNP(file_name):
    # Load dataset
    ds = xr.open_dataset(file_name)
    # Create mask based on Basin_mask and Harvest Area threshold
    mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, 0)
    
    # Calculate masked variables
    Total_Irri = ds["Total_irrigation_amount"] * mask
    Sus_Irri   = ds["Sus_irrigation_amount"] * mask
    N_runoff   = ds["N_Runoff"] * mask
    P_runoff   = ds["P_Runoff"] * mask
    Crit_N     = ds["Crit_N_Runoff"] * mask
    Crit_P     = ds["Crit_P_Runoff"] * mask
    
    return Total_Irri, Sus_Irri, N_runoff, P_runoff, Crit_N, Crit_P

for basin in Studyarea:
    print(f"Processing Basin: {basin}")
    
    # Initialize basin-level sums as 0
    basin_W = 0; basin_W_crit = 0
    basin_N = 0; basin_N_crit = 0
    basin_P = 0; basin_P_crit = 0
    
    # Track if we found any valid data for this basin
    found_data = False

    for crop in Croptypes:
        # Construct path for each crop file
        # Note: adjust the filename pattern if it differs from {basin}_{crop}_summary.nc
        file_name = os.path.join(model_output_dir, f"{basin}_{crop}_summary.nc")
        
        if not os.path.exists(file_name):
            continue 
        
        # Extract and sum
        W, Wc, N, P, Nc, Pc = GetCropWNP(file_name)
        
        basin_W += W.fillna(0)
        basin_W_crit += Wc.fillna(0)
        basin_N += N.fillna(0)
        basin_N_crit += Nc.fillna(0)
        basin_P += P.fillna(0)
        basin_P_crit += Pc.fillna(0)
        found_data = True

    if not found_data:
        print(f"No crop data found for {basin}, skipping...")
        continue

    # --- Load Low Runoff Mask ---
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]

    # --- CLASSIFICATION LOGIC ---
    # We use boolean logic to build the codes (Safe=1, Exceeded=2)
    # W_digit: 100 if Safe, 200 if Exceeded
    # N_digit: 10 if Safe, 20 if Exceeded
    # P_digit: 1 if Safe, 2 if Exceeded
    
    W_exceed = xr.where(basin_W > basin_W_crit, 200, 100)
    N_exceed = xr.where(basin_N > basin_N_crit, 20, 10)
    P_exceed = xr.where(basin_P > basin_P_crit, 2, 1)
    
    # Combine to get 111, 112, 121, etc.
    pixel_value = W_exceed + N_exceed + P_exceed
    
    # Priority 1: If no crop activity (Sum of runoff is 0), set value to 0
    # Priority 2: If low_runoff == 1, set value to 1
    final_mask = pixel_value.where((basin_N + basin_P + basin_W) > 0, 0)
    final_mask = xr.where(low_runoff == 1, 1, final_mask)

    # --- Save Output ---
    # Ensure dimensions are named correctly for GeoTIFF (usually 'lat'/'lon' to 'y'/'x')
    final_mask_rio = final_mask.rename({'lon': 'x', 'lat': 'y'})
    
    # Assign the coordinate reference system (WGS84)
    final_mask_rio.rio.write_crs("epsg:4326", inplace=True)

    final_mask_rio.to_netcdf(os.path.join(output_dir, f"{basin}_boundary_check.nc"))

    print(f"Saved NetCDF for {basin}")      