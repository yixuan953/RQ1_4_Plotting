# This script is used plot how do N runoff, N uptake, N in grain, Crop production change with different scenarios.

import os
import numpy as np
import xarray as xr
import pandas as pd

Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
InputCrops = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

Irrigation_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline"
Rainfed_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline_rainfed"

Irrigated_fert_red_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert"
Rainfed_fert_red_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed"

Scenarios = ["Red_02", "Red_04", "Red_06", "Red_08", "Red_10", "Red_12", "Red_14", "Red_prop", "Inc_02", "Inc_04", "Inc_06", "Inc_08", "Inc_10", "Inc_12", "Inc_14", "Inc_prop"]

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
fert_red_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/4_Fertilization_Red/4_1_Reduced_Fert"

csv_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/4_Fertilization_Red/4_2_Summary"
HA_Threshold = 2500

def ReadModelOutput(file_path):
    if not os.path.exists(file_path):
        return None
    ds = xr.open_dataset(file_path)
    # Return the variables slice-only; don't mask here yet
    return (ds["N_Runoff"].sel(year=slice(2010, 2019)), 
            ds["N_Uptake"].sel(year=slice(2010, 2019)), 
            ds["N_grain"].sel(year=slice(2010, 2019)), 
            ds["Yield"].sel(year=slice(2010, 2019)))

def ReadFertilizer(file_path, mask_with_HA, mask_not_low_runoff, cropname):
    if not os.path.exists(file_path):
        print(f"File {file_path} does not exist. Skipping...")
        return None
    ds = xr.open_dataset(file_path)
    Inorg = ds[f"{cropname}_Inorg_N_application_rate"].sel(year=slice(2010, 2019))
    Urea = ds[f"{cropname}_Urea_N_application_rate"].sel(year=slice(2010, 2019))
    Manure = ds[f"{cropname}_Manure_N_application_rate"].sel(year=slice(2010, 2019))
    Fertilizer = Inorg + Urea + Manure
    Fertilizer = Fertilizer * mask_with_HA * mask_not_low_runoff
    return Fertilizer

def BasinStatistics(data_var, HA_after_masking, mask_not_low_runoff):
    # Align both masks to the data_var grid
    mask_lr = mask_not_low_runoff.reindex_like(data_var, method="nearest")
    ha_aligned = HA_after_masking.reindex_like(data_var, method="nearest")
    
    # Create the final weight: Area where it's NOT low runoff
    weights = ha_aligned * mask_lr
    
    Total_amount = (data_var * weights).sum(dim=['lat', 'lon'], skipna=True)
    Total_harvested_area = weights.sum(dim=['lat', 'lon'], skipna=True)
    
    return xr.where(Total_harvested_area > 0, Total_amount/Total_harvested_area, 0.0)

# Initialize a list to store results

for basin in Studyareas:
    # Load Masks once per basin
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        # Assuming the coord names are 'lat' and 'lon'
        mask_not_low_runoff = xr.where(ds_lr["Low_Runoff"].isnull(), 1, np.nan)

    results_list = []

    for crop in InputCrops:

        if crop == "winterwheat": cropname1 = cropname2 = "Wheat"
        elif crop == "maize": cropname1 = cropname2 = "Maize"
        elif crop == "mainrice": cropname1 = cropname2 = "Rice"
        elif crop == "secondrice": cropname1 = "Rice"; cropname2 = "Secondrice"
        elif crop == "soybean": cropname1 = cropname2 = "Soybean"

        summary_path = os.path.join(summary_dir, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(summary_path):
            continue

        print(f"Processing {basin} - {crop}...")
        
        ds_summary = xr.open_dataset(summary_path)
        total_HA = ds_summary["Total_HA"].where(ds_summary["Total_HA"] > HA_Threshold, 0)
        mask_with_HA = xr.where(total_HA > HA_Threshold, 1, np.nan)

        Irrigated_HA = ds_summary["Irrigated_HA"].where(ds_summary["Irrigated_HA"] > HA_Threshold, 0)
        Rainfed_HA = ds_summary["Rainfed_HA"].where(ds_summary["Rainfed_HA"] > HA_Threshold, 0)
        
        # Combine Scenarios with Baseline for a single loop
        all_scenarios = ["Baseline"] + Scenarios

        for sc in all_scenarios:
            # Set up paths based on whether it is Baseline or a specific reduction scenario
            if sc == "Baseline":
                path_map = {
                    "Irrigated": (os.path.join(Irrigation_baseline_dir, f"{basin}_{crop}_annual.nc"), 
                                os.path.join(data_dir, "2_StudyArea", basin, "Fertilization", f"{basin}_{cropname1}_Fert_2005-2020_FixRate.nc"), Irrigated_HA),
                    "Rainfed": (os.path.join(Rainfed_baseline_dir, f"{basin}_{crop}_annual.nc"),
                               os.path.join(data_dir, "2_StudyArea", basin, "Fertilization", f"{basin}_{cropname1}_Fert_2005-2020_FixRate.nc"), Rainfed_HA)
                }
            else:
                path_map = {
                    "Irrigated": (os.path.join(Irrigated_fert_red_dir, sc, f"{basin}_{crop}_annual.nc"),
                                 os.path.join(fert_red_dir, "Irrigated", sc, f"{basin}_{cropname2}_Fert_2005-2020_FixRate.nc"), Irrigated_HA),
                    "Rainfed": (os.path.join(Rainfed_fert_red_dir, sc, f"{basin}_{crop}_annual.nc"),
                               os.path.join(fert_red_dir, "Rainfed", sc, f"{basin}_{cropname2}_Fert_2005-2020_FixRate.nc"), Rainfed_HA)
                }

            for system_type, (model_path, fert_path, Harvested_Area) in path_map.items():
                # Read Data
                model_data = ReadModelOutput(model_path)
                fert_data = ReadFertilizer(fert_path, mask_with_HA, mask_not_low_runoff, cropname1)

                if model_data is None or fert_data is None:
                    continue

                n_run, n_upt, n_grain, c_prod = model_data
                
                # Calculate Basin Averages (Returns a time-series per year)
                avg_run = BasinStatistics(n_run, Harvested_Area, mask_not_low_runoff)
                avg_upt = BasinStatistics(n_upt, Harvested_Area, mask_not_low_runoff)
                avg_grain = BasinStatistics(n_grain, Harvested_Area, mask_not_low_runoff)
                avg_prod = BasinStatistics(c_prod, Harvested_Area, mask_not_low_runoff)
                avg_fert = BasinStatistics(fert_data, Harvested_Area, mask_not_low_runoff)

                # Loop through years 2010-2019 to extract values
                for year in range(2010, 2020):
                    results_list.append({
                        "Scenario": sc,
                        "Year": year,
                        "Crop": cropname2,
                        "System": system_type,
                        "Fert": float(avg_fert.sel(year=year)),
                        "N_runoff": float(avg_run.sel(year=year)),
                        "N_uptake": float(avg_upt.sel(year=year)),
                        "N_grain": float(avg_grain.sel(year=year)),
                        "Yield": float(avg_prod.sel(year=year))
                    })

    # Finalize and Save
    df = pd.DataFrame(results_list)
    os.makedirs(csv_output_dir, exist_ok=True)
    df.to_csv(os.path.join(csv_output_dir, f"{basin}_basin_Avg_Summary.csv"), index=False)
    print(f"{basin}_basin_Avg_Summary.csv saved.")
        