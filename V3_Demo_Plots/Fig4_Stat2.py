# This script is used to calculate how much fertilization input has been reduced compared to the baseline scenario each crop type (Rainfed, Irrigated, and Total):
# And save the output to .csv

# 1 - Total irrigated area [ha]
# 2 - Total rainfed area [ha]
# 3 - Total harvested area [ha]

# 4 - Current_N_Inorg  [ktons]
# 5 - Current_N_Manure [ktons]
# 6 - Current_P_Inorg  [ktons]
# 7 - Current_P_Manure [ktons]

# 8 - Reduced_N_Inorg_Rainfed [ktons]
# 9 - Reduced_N_Manure_Rainfed [ktons]
# 10 - Reduced_P_Inorg_Rainfed [ktons]
# 11 - Reduced_P_Manure_Rainfed [ktons]

# 12 - Reduced_N_Inorg_Irrigated [ktons]
# 13 - Reduced_N_Manure_Irrigated [ktons]
# 14 - Reduced_P_Inorg_Irrigated [ktons]
# 15 - Reduced_P_Manure_Irrigated [ktons]

# 16 - Current_N_Inorg_Fert_Rate  [kg/ha]
# 17 - Current_N_Manure_Fert_Rate  [kg/ha]
# 18 - Current_P_Inorg_Fert_Rate  [kg/ha]
# 19 - Current_P_Manure_Fert_Rate  [kg/ha]

# 20 - Reduced_N_Inorg_Rainfed_Fert_Rate  [kg/ha]
# 21 - Reduced_N_Manure_Rainfed_Fert_Rate  [kg/ha]
# 22 - Reduced_P_Inorg_Rainfed_Fert_Rate  [kg/ha]
# 23 - Reduced_P_Manure_Rainfed_Fert_Rate  [kg/ha]

# 24 - Reduced_N_Inorg_Irrigated_Fert_Rate  [kg/ha]
# 25 - Reduced_N_Manure_Irrigated_Fert_Rate  [kg/ha]
# 26 - Reduced_P_Inorg_Irrigated_Fert_Rate  [kg/ha]
# 27 - Reduced_P_Manure_Irrigated_Fert_Rate  [kg/ha]

import os
import numpy as np
import xarray as xr
import pandas as pd

# --- Paths ---
model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
baseline_fert_input_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
Irrigated_fert_input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/4_Fertilization_Red/4_1_Reduced_Fert/Irrigated/Red_prop"
Rainfed_fert_input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/4_Fertilization_Red/4_1_Reduced_Fert/Rainfed/Red_prop"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/4_Red_Fert"
os.makedirs(output_dir, exist_ok=True)

STUDY_AREAS = ["LaPlata", "Rhine", "Indus", "Yangtze"]
# Exact names for file pathing
CROP_GROUPS = ["Wheat", "Maize", "Rice", "Soybean"]

def ReadHarvestedArea(file_name, low_runoff_path):
    if not os.path.exists(file_name): return None
    with xr.open_dataset(file_name) as ds:
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
        with xr.open_dataset(low_runoff_path) as ds_lr:
            low_runoff = ds_lr["Low_Runoff"].reindex_like(mask, method='nearest')
        mask = mask.where(low_runoff != 1, np.nan)
        
        return {
            "total": ds["Total_HA"] * mask,
            "rain": ds["Rainfed_HA"] * mask,
            "irri": ds["Irrigated_HA"] * mask
        }, mask

def GetFertMass(nc_path, crop_fn, area_mask):
    """Calculates total mass in ktons (kg / 1e6) using the rate file."""
    if not os.path.exists(nc_path): return 0, 0, 0, 0
    with xr.open_dataset(nc_path) as ds:
        # Select 2015 and align to mask
        ds_2015 = ds.sel(year=2015).reindex_like(area_mask, method='nearest')
        
        # Mass [kg] = Rate [kg/ha] * Area [ha]
        if crop_fn == "Secondrice":
            crop_fn2 = "Rice"
        else:
            crop_fn2 = crop_fn
        m_ni = (ds_2015[f"{crop_fn2}_Total_inorg_N_application_rate"] * area_mask).sum().values * 1e-6
        m_nm = (ds_2015[f"{crop_fn2}_Manure_N_application_rate"] * area_mask).sum().values * 1e-6 
        m_pi = (ds_2015[f"{crop_fn2}_P_application_rate"] * area_mask).sum().values * 1e-6
        m_pm = (ds_2015[f"{crop_fn2}_Manure_P_application_rate"] * area_mask).sum().values * 1e-6
        
    return float(m_ni), float(m_nm), float(m_pi), float(m_pm)

# --- Main Loop ---
for basin in STUDY_AREAS:
    low_runoff = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    basin_results = []
    
    for crop in CROP_GROUPS:
        # Determine files to process for this group
        sub_crops = ["Rice", "Secondrice"] if (basin == "Yangtze" and crop == "Rice") else [crop]
        
        # Accumulators for the group (e.g. combined Rice + Secondrice)
        group_ha_t, group_ha_r, group_ha_i = 0, 0, 0
        group_base_ni, group_base_nm, group_base_pi, group_base_pm = 0, 0, 0, 0
        group_sc_ni, group_sc_nm, group_sc_pi, group_sc_pm = 0, 0, 0, 0

        for sc in sub_crops:
            mapping = {
                "Wheat": "winterwheat",
                "Maize": "maize",
                "Rice": "mainrice",
                "Secondrice": "secondrice",
                "Soybean": "soybean"
            }
            summary_fn = mapping.get(sc)
                
            nc_summary = os.path.join(model_summary_dir, f"{basin}_{summary_fn}_summary.nc")
            if not os.path.exists(nc_summary):
                print(f"Warning: Summary file not found: {nc_summary}")
                continue
            
            area_data = ReadHarvestedArea(nc_summary, low_runoff)
            if area_data is None: continue
            areas, mask = area_data
            
            # Paths to fertilization rate files
            base_f = os.path.join(baseline_fert_input_dir, f"{basin}/Fertilization/{basin}_{sc}_Fert_2005-2020_FixRate.nc")
            irri_f = os.path.join(Irrigated_fert_input_dir, f"{basin}_{sc}_Fert_2005-2020_FixRate.nc")
            rain_f = os.path.join(Rainfed_fert_input_dir, f"{basin}_{sc}_Fert_2005-2020_FixRate.nc")

            # 1. Sum Areas
            group_ha_t += float(areas['total'].sum())
            group_ha_r += float(areas['rain'].sum())
            group_ha_i += float(areas['irri'].sum())

            # 2. Get Baseline Mass (Total Area)
            bni, bnm, bpi, bpm = GetFertMass(base_f, sc, areas['total'])
            group_base_ni += bni; group_base_nm += bnm; group_base_pi += bpi; group_base_pm += bpm

            # 3. Get Post-Reduction Mass
            # We need the Scenario rates on their respective areas
            sni_r, snm_r, spi_r, spm_r = GetFertMass(rain_f, sc, areas['rain'])
            sni_i, snm_i, spi_i, spm_i = GetFertMass(irri_f, sc, areas['irri'])
            
            group_sc_ni += (sni_r + sni_i)
            group_sc_nm += (snm_r + snm_i)
            group_sc_pi += (spi_r + spi_i)
            group_sc_pm += (spm_r + spm_i)

        if group_ha_t == 0: continue

        # 4. Final Aggregation (Current vs Post-Reduction)
        res = {
            "Crop_Group": crop,
            "Total_HA": group_ha_t,
            "Rainfed_HA": group_ha_r,
            "Irrigated_HA": group_ha_i,

            # Current & Post-Redudction Totals [ktons]
            "Current_N_Inorg_kt": group_base_ni,
            "Current_N_Manure_kt": group_base_nm,
            "PostRed_N_Inorg_kt": group_sc_ni,
            "PostRed_N_Manure_kt": group_sc_nm,

            "Current_P_Manure_kt": group_base_pm,
            "Current_P_Inorg_kt": group_base_pi,
            "PostRed_P_Inorg_kt": group_sc_pi,
            "PostRed_P_Manure_kt": group_sc_pm,

            # Average Rates [kg/ha] (Combined Mass / Combined Area)
            "Current_Rate_N_Inorg": (group_base_ni * 1e6 / group_ha_t),
            "Current_Rate_N_Manure": (group_base_nm * 1e6 / group_ha_t),
            "PostRed_Rate_N_Inorg": (group_sc_ni * 1e6 / group_ha_t),
            "PostRed_Rate_N_Manure": (group_sc_nm * 1e6 / group_ha_t),

            "Current_Rate_P_Inorg": (group_base_pi * 1e6 / group_ha_t),
            "Current_Rate_P_Manure": (group_base_pm * 1e6 / group_ha_t),
            "PostRed_Rate_P_Inorg": (group_sc_pi * 1e6 / group_ha_t),
            "PostRed_Rate_P_Manure": (group_sc_pm * 1e6 / group_ha_t)
        }
        basin_results.append(res)

    if basin_results:
        pd.DataFrame(basin_results).to_csv(os.path.join(output_dir, f"{basin}_Fert_Final_Summary.csv"), index=False)

print("Calculations finished.")

        
