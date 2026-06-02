# This script is used to calculate how many ha of harvested area and cropland area has exceeded each boundary (irrigation, N runoff, P runoff) for each crop and each basin

import pandas as pd
import numpy as np
import xarray as xr
import os

# model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
# output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Statistics/3_Exceedance"

model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Statistics/5_Implications"

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

os.makedirs(output_dir, exist_ok=True)

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

def GetCropWNP(file_name, low_runoff_path, PA_path):

    with xr.open_dataset(file_name) as ds:
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)

    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]
    
    mask = mask.where(low_runoff != 1, np.nan)

    with xr.open_dataset(PA_path) as ds_pa:
        Total_PA = ds_pa["Physical_Area"] * mask
        Irri_PA = ds_pa["Physical_Area"] * ds_pa["Irrigated_Proportion"] * mask
    
    # Extracting variables
    Total_HA = ds["Total_HA"] * mask
    Irri_HA = ds["Irrigated_HA"] * mask
    T_Irr = ds["Total_irrigation_amount"] * mask
    S_Irr = ds["Sus_irrigation_amount"] * mask
    N = ds["N_Runoff"] * mask
    CN = ds["Crit_N_Runoff"] * mask
    P = ds["P_Runoff"] * mask
    CP = ds["Crit_P_Runoff"] * mask

    return Total_HA, Irri_HA, T_Irr, S_Irr, N, P, CN, CP, Total_PA, Irri_PA

# Loop through each basin to create individual CSVs
for basin in Studyarea:

    basin_results_ha = []
    basin_results_pa = []
    
    for crop in Croptypes:
        file_name = os.path.join(model_summary_dir, f"{basin}_{crop}_summary.nc")
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")

        if crop == "winterwheat":
            cropname = "WHEA"
        elif crop == "maize":
            cropname = "MAIZ"
        elif crop == "mainrice":
            cropname = "RICE"
        elif crop == "secondrice":
            cropname = "RICE"
        elif crop == "soybean":
            cropname = "SOYB"
        PA_path = os.path.join(data_dir, "2_StudyArea", basin, f"PhysicalArea/{cropname}_PA_05d.nc")
        
        if os.path.exists(file_name) and os.path.exists(low_runoff_path) and os.path.exists(PA_path):
            HA, Irri_HA, T_Irr, S_Irr, N, P, CN, CP, Total_PA, Irri_PA = GetCropWNP(file_name, low_runoff_path, PA_path)
            
            # Boolean masks
            w_exc = (T_Irr > S_Irr)
            n_exc = (N > CN)
            p_exc = (P > CP)

            total_ha_sum = float(HA.sum())
            total_pa_sum = float(Total_PA.sum())
            
            if total_ha_sum > 0:
                row_ha = {
                    "Crop type": crop,
                    "Irri": float(Irri_HA.where(w_exc).sum()) / total_ha_sum,
                    "N": float(HA.where(n_exc).sum()) / total_ha_sum,
                    "P": float(HA.where(p_exc).sum()) / total_ha_sum,
                    "Irri & N": float(Irri_HA.where(w_exc & n_exc).sum()) / total_ha_sum,
                    "Irri & P": float(Irri_HA.where(w_exc & p_exc).sum()) / total_ha_sum,
                    "N & P": float(HA.where(n_exc & p_exc).sum()) / total_ha_sum,
                    "All three": float(Irri_HA.where(w_exc & n_exc & p_exc).sum()) / total_ha_sum
                }
                basin_results_ha.append(row_ha)

                row_pa = {
                    "Crop type": crop,
                    "Irri": float(Irri_PA.where(w_exc).sum()) / total_pa_sum,
                    "N": float(Total_PA.where(n_exc).sum()) / total_pa_sum,
                    "P": float(Total_PA.where(p_exc).sum()) / total_pa_sum,
                    "Irri & N": float(Irri_PA.where(w_exc & n_exc).sum()) / total_pa_sum,
                    "Irri & P": float(Irri_PA.where(w_exc & p_exc).sum()) / total_pa_sum,
                    "N & P": float(Total_PA.where(n_exc & p_exc).sum()) / total_pa_sum,
                    "All three": float(Irri_PA.where(w_exc & n_exc & p_exc).sum()) / total_pa_sum
                }
                basin_results_pa.append(row_pa)

    # Process the results for this specific basin
    if basin_results_ha:
        df_basin_ha = pd.DataFrame(basin_results_ha)
        output_file_ha = os.path.join(output_dir, f"{basin}_boundaries_ExceedanceHA.csv")
        df_basin_ha.to_csv(output_file_ha, index=False)
        print(f"Saved: {output_file_ha}")

    if basin_results_pa:
        df_basin_pa = pd.DataFrame(basin_results_pa)
        output_file_pa = os.path.join(output_dir, f"{basin}_boundaries_ExceedancePA.csv")
        df_basin_pa.to_csv(output_file_pa, index=False)
        print(f"Saved: {output_file_pa}")