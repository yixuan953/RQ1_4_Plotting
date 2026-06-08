# This script is used to calculate the average, median, and 25th and 75th yield in the basin for years (2010 - 2019)

import pandas as pd
import numpy as np
import os
import xarray as xr

# Read the summary file
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/SM/1_Val_CropYield"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

os.makedirs(output_dir, exist_ok=True)

studyareas = ["Indus", "Rhine", "LaPlata", "Yangtze"]
croptypes = ["mainrice", "secondrice", "maize", "soybean", "winterwheat"]

def weighted_stats(values, weights):
    """Return weighted mean, median, p25, p75"""
    values = np.array(values)
    weights = np.array(weights)

    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[mask]
    weights = weights[mask]

    if len(values) == 0:
        return np.nan, np.nan, np.nan, np.nan

    # weighted mean
    w_mean = np.sum(values * weights) / np.sum(weights)

    # weighted quantiles
    sorter = np.argsort(values)
    values_sorted = values[sorter]
    weights_sorted = weights[sorter]
    cum_w = np.cumsum(weights_sorted) / np.sum(weights_sorted)

    def weighted_quantile(q):
        return np.interp(q, cum_w, values_sorted)

    w_median = weighted_quantile(0.5)
    p25 = weighted_quantile(0.25)
    p75 = weighted_quantile(0.75)

    return w_mean, w_median, p25, p75


results = []

for basin in studyareas:
    for crop in croptypes:

        input_file = f"{input_dir}/{basin}_{crop}_summary.nc"
        if not os.path.exists(input_file):
            print(f"!!! {input_file} missing, skipping...")
            continue

        with xr.open_dataset(input_file) as ds:

            # --- variables ---
            rain_yield = ds["Avg_Yield_Rainfed"].values
            irr_yield = ds["Avg_Yield_Irrigated"].values
            rain_area = ds["Avg_Rainfed_Area"].values
            irr_area = ds["Avg_Irrigated_Area"].values
            basin_mask = ds["Basin_mask"].values  # important: keep values outside context

        # --- low runoff mask ---
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")

        with xr.open_dataset(low_runoff_path) as ds_lr:
            low_runoff = ds_lr["Low_Runoff"].values

        # --- combine masks ---
        valid_mask = (basin_mask == 1) & (low_runoff != 1)

        # flatten
        rain_yield = rain_yield[valid_mask]
        irr_yield = irr_yield[valid_mask]
        rain_area = rain_area[valid_mask]
        irr_area = irr_area[valid_mask]

        # --- compute stats ---
        rain_stats = weighted_stats(rain_yield, rain_area)
        irr_stats = weighted_stats(irr_yield, irr_area)

        results.append({
            "Basin": basin,
            "Crop": crop,

            "Rainfed_mean": rain_stats[0],
            "Rainfed_median": rain_stats[1],
            "Rainfed_p25": rain_stats[2],
            "Rainfed_p75": rain_stats[3],

            "Irrigated_mean": irr_stats[0],
            "Irrigated_median": irr_stats[1],
            "Irrigated_p25": irr_stats[2],
            "Irrigated_p75": irr_stats[3],
        })

# --- save ---
stats_df = pd.DataFrame(results)

output_file = os.path.join(output_dir, "Simulated_Crop_Yield.csv")
stats_df.to_csv(output_file, index=False)

print("Done:", output_file)