# This script is used to calculate
# 1 - Total exceedance of N boundaries by crop type [ktons]
# 2 - Total exceedance of P boundaries by crop type  [ktons]
# 3 - Basin average exceedance of N boundaries by crop type [kg/ha]
# 4 - Basin average exceedance of P boundaries by crop type [kg/ha]

# Save the output to .csv file

import pandas as pd
import numpy as np
import xarray as xr
import os

model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/3_Exceedance"

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
row_labels = ["Total_exceedance_N", "Total_exceedance_P", "Avg_N_exceedance", "Avg_P_exceedance"]

def GetCropWNP(file_name, low_runoff_path):
    with xr.open_dataset(file_name) as ds:
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)

    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]
        
    mask = mask.where(low_runoff != 1, np.nan)
        
    # 1. Hectares
    t_ha = float((ds["Total_HA"] * mask).sum())

    # 2. Critical Loads
    CN = (ds["Crit_N_Runoff"] * mask).sum().values
    CP = (ds["Crit_P_Runoff"] * mask).sum().values 
    # 3. Total Runoff (Current)
    total_N = float(((ds["N_Runoff_Irrigated"] * ds["Irrigated_HA"] * mask + ds["N_Runoff_Rainfed"] * ds["Rainfed_HA"]) * mask).sum())
    total_P = float(((ds["P_Runoff_Irrigated"] * ds["Irrigated_HA"] * mask + ds["P_Runoff_Rainfed"] * ds["Rainfed_HA"]) * mask).sum())

    # 4. Exceedance
    Total_N_exceedance = total_N - CN
    Total_P_exceedance = total_P - CP 

    # 5. Basin Average [kg/ha]
    Avg_N_exceedance = Total_N_exceedance / t_ha if t_ha > 0 else 0
    Avg_P_exceedance = Total_P_exceedance / t_ha if t_ha > 0 else 0

    # Convert Totals to ktons (1e-6)
    Total_N_exceedance_kt = Total_N_exceedance * 1e-6
    Total_P_exceedance_kt = Total_P_exceedance * 1e-6

    stats = [Total_N_exceedance_kt, Total_P_exceedance_kt, Avg_N_exceedance, Avg_P_exceedance]
    areas = [t_ha]
    
    return stats, areas

# --- Main Loop ---
for basin in Studyarea:
    basin_results = {}
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
            basin_results[display_name] = [np.nan] * 4
            continue

        if len(all_stats) > 1:
            # Indices 0, 1 are TOTALS [ktons] -> Sum them
            # Indices 2, 3 are AVERAGES [kg/ha] -> Area-weighted average
            t_n_ex = sum(item[0] for item in all_stats)
            t_p_ex = sum(item[1] for item in all_stats)
            
            total_area = sum(a[0] for a in all_areas)
            
            if total_area > 0:
                # Weighted avg = Sum(Avg * Area) / Total Area
                avg_n_ex = sum(all_stats[j][2] * all_areas[j][0] for j in range(len(all_stats))) / total_area
                avg_p_ex = sum(all_stats[j][3] * all_areas[j][0] for j in range(len(all_stats))) / total_area
            else:
                avg_n_ex, avg_p_ex = 0, 0
                
            basin_results[display_name] = [t_n_ex, t_p_ex, avg_n_ex, avg_p_ex]
        else:
            basin_results[display_name] = all_stats[0]

    # Save CSV - Using float_format to keep integers in the output
    df = pd.DataFrame(basin_results, index=row_labels)
    output_path = os.path.join(output_dir, f"{basin}_exceedance_summary.csv")
    df.to_csv(output_path, float_format="%.2f")
    print(f"Saved: {output_path}")
            