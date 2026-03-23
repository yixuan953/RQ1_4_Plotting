# This script is used to calculate
# 1 - Total crop production of rainfed & irrigated field [kg]
# 2 - Total cropland N runoff to surface runoff of rainfed & irrigated field [kg]
# 3 - Total cropland P runoff to surface runoff of rainfed & irrigated field [kg]

# 5 - Basin average crop yield of rainfed & irrigated field [kg/ha]
# 6 - Basin average cropland N runoff to surface runoff of rainfed & irrigated field [kg/ha]
# 7 - Basin average cropland P runoff to surface runoff of rainfed & irrigated field [kg/ha]

# Save the output to .csv file

import pandas as pd
import numpy as np
import xarray as xr
import os

model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/2_Simulated_Yield_Runoff/1_Current"

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

row_labels = [
    "Total_Production_Rainfed", "Total_Production_Irrigated", 
    "Total_N_runoff_Rainfed", "Total_N_runoff_Irrigated", 
    "Total_P_runoff_Rainfed", "Total_P_runoff_Irrigated",
    "Avg_yield_rainfed", "Avg_yield_irrigated", "Avg_yield_total",
    "Avg_N_runoff_rainfed", "Avg_N_runoff_irrigated", "Avg_N_runoff_total",
    "Avg_P_runoff_rainfed", "Avg_P_runoff_irrigated", "Avg_P_runoff_total"
]

def GetCropWNP(file_name, low_runoff_path):
    with xr.open_dataset(file_name) as ds:
        mask = xr.where(ds["Total_HA"] > 2500, 1, np.nan)

    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]
        
    mask = mask.where(low_runoff != 1, np.nan)
        
    # 1. Hectares (Basin Totals)
    t_ha = float((ds["Total_HA"] * mask).sum())
    t_ha_r = float((ds["Rainfed_HA"] * mask).sum())
    t_ha_i = float((ds["Irrigated_HA"] * mask).sum())

    # 2. Production [ktons] (Sum the grid first, then convert to float)
    total_Prod_Irri = float((ds["Avg_Yield_Irrigated"] * ds["Irrigated_HA"] * mask).sum())
    total_Prod_Rain = float((ds["Avg_Yield_Rainfed"] * ds["Rainfed_HA"] * mask).sum())

    # 3. Runoff [ktons] (Sum the grid first)
    total_N_Irri = float((ds["N_Runoff_Irrigated"] * ds["Irrigated_HA"]  * mask).sum())
    total_N_Rain = float((ds["N_Runoff_Rainfed"] * ds["Rainfed_HA"] * mask).sum())
    
    total_P_Irri = float((ds["P_Runoff_Irrigated"] * ds["Irrigated_HA"]  * mask).sum())
    total_P_Rain = float((ds["P_Runoff_Rainfed"] * ds["Rainfed_HA"]  * mask).sum())

    # 4. Basin Averages [kg/ha] 
    # Logic: (ktons * 1e6) / ha
    avg_y_r = (total_Prod_Rain / t_ha_r) if t_ha_r > 0 else 0
    avg_y_i = (total_Prod_Irri / t_ha_i) if t_ha_i > 0 else 0
    avg_y_t = ((total_Prod_Rain + total_Prod_Irri) / t_ha) if t_ha > 0 else 0

    avg_N_r = (total_N_Rain / t_ha_r) if t_ha_r > 0 else 0
    avg_N_i = (total_N_Irri / t_ha_i) if t_ha_i > 0 else 0
    avg_N_t = ((total_N_Rain + total_N_Irri) / t_ha) if t_ha > 0 else 0

    avg_P_r = (total_P_Rain / t_ha_r) if t_ha_r > 0 else 0
    avg_P_i = (total_P_Irri / t_ha_i) if t_ha_i > 0 else 0
    avg_P_t = ((total_P_Rain + total_P_Irri) / t_ha) if t_ha > 0 else 0

    # We return the 15 stats, AND the 3 area values to use for combining crops later
    stats = [total_Prod_Rain, total_Prod_Irri, total_N_Rain, total_N_Irri, total_P_Rain, total_P_Irri,
             avg_y_r, avg_y_i, avg_y_t, avg_N_r, avg_N_i, avg_N_t, avg_P_r, avg_P_i, avg_P_t]
    areas = [t_ha_r, t_ha_i, t_ha]
    
    return stats, areas

# --- Main Loop ---

for basin in Studyarea:
    basin_results = {}
    # Crop logic: Yangtze has two rice crops, others have one
    crop_map = {
        "Wheat": ["winterwheat"], 
        "Maize": ["maize"], 
        "Rice": ["mainrice", "secondrice"] if basin == "Yangtze" else ["mainrice"], 
        "Soybean": ["soybean"]
    }

    for display_name, sub_crops in crop_map.items():
        all_stats = []
        all_areas = []
        
        for crop in sub_crops:
            nc_summary_file = os.path.join(model_summary_dir, f"{basin}_{crop}_summary.nc")
            low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")

            if os.path.exists(nc_summary_file):
                s, a = GetCropWNP(nc_summary_file, low_runoff_path)
                all_stats.append(s)
                all_areas.append(a)

        if not all_stats:
            basin_results[display_name] = [np.nan] * 15
            continue

        # COMBINING CROPS (e.g., Mainrice + Secondrice)
        if len(all_stats) > 1:
            combined = []
            # Indices 0-5 are TOTALS (Sum them directly)
            for i in range(6):
                combined.append(sum(item[i] for item in all_stats))
            
            # Indices 6-14 are AVERAGES (Recalculate using Area weighting)
            # Area mapping: indices 6,9,12 use Rainfed_HA [0]; 7,10,13 use Irri_HA [1]; 8,11,14 use Total_HA [2]
            area_weight_idx = [0, 1, 2] * 3 
            for i in range(6, 15):
                w_idx = area_weight_idx[i-6]
                total_weighted_val = sum(all_stats[j][i] * all_areas[j][w_idx] for j in range(len(all_stats)))
                total_area = sum(all_areas[j][w_idx] for j in range(len(all_areas)))
                combined.append(total_weighted_val / total_area if total_area > 0 else 0)
            
            basin_results[display_name] = combined
        else:
            basin_results[display_name] = all_stats[0]

    # Save CSV
    df = pd.DataFrame(basin_results, index=row_labels)
    df.to_csv(os.path.join(output_dir, f"{basin}_summary_stats.csv"))
            