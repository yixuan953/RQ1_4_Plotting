# This script is used to output the annual basin average p pool and P runoff in .csv file

import os
import xarray as xr
import numpy as np
import pandas as pd

# Configuration
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
CropTypes = ["winterwheat", "mainrice", "secondrice", "maize", "soybean"]
Sens_Scenarios = ["Red_prop", "Red_02", "Red_04", "Red_06", "Red_08", "Red_10", "Red_12", "Red_14"] 

p_mol = 30.97 # 1 mol P = 30.97 g

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig5/5c"
os.makedirs(output_dir, exist_ok=True)

def get_metrics(csv_path, area_df, bd_df):
    if not os.path.exists(csv_path):
        return None, None
    
    df = pd.read_csv(csv_path, skipinitialspace=True)
    df = df.rename(columns={'Lat': 'lat', 'Lon': 'lon'})
    df = df[(df["Year"] >= 2009) & (df["Year"] <= 2019)]
    
    # Merge datasets
    m = df.merge(area_df, on=["lat", "lon"]).merge(bd_df, on=["lat", "lon"])
    
    # Calculate weighted components
    m["pixel_mass"] = m["HA"] * m["bulk_density"]
    m["weighted_P"] = (m["LabileP"] + m["StableP"]) * m["pixel_mass"] * p_mol
    m["weighted_P_runoff"] = (m["P_surf"] + m["P_sub"]) * m["HA"]
    
    grouped = m.groupby("Year")
    
    # We return the mass as well to correctly re-weight during total aggregation
    total_mass = grouped["pixel_mass"].sum()
    total_HA = grouped["HA"].sum()
    p_pool = grouped["weighted_P"].sum() / total_mass
    p_runoff = grouped["weighted_P_runoff"].sum() / total_HA
    
    return p_pool, p_runoff, total_mass

for basin in Studyareas:
    range_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "range.nc")
    
    with xr.open_dataset(range_nc) as ds_range:
        template = ds_range["mask"]

    low_runoff_path = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan) 

    for crop in CropTypes:
        # 1 - Paths for Harvested area
        sc1_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline", f"{basin}_{crop}_summary.nc")
        if not os.path.exists(sc1_nc):
            continue
            
        with xr.open_dataset(sc1_nc) as ds:
            mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan) * template * mask_not_low_runoff
            irri_ha = (ds["Irrigated_HA"] * mask).to_dataframe(name="HA").reset_index().dropna()
            rain_ha = (ds["Rainfed_HA"] * mask).to_dataframe(name="HA").reset_index().dropna()

        # 2 - Path for soil bulk density
        soil_prop_file = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "Mask", f"{basin}_{crop}_mask.nc")
        if not os.path.exists(soil_prop_file):
            continue
        with xr.open_dataset(soil_prop_file) as soil_ds:
            bd_df = (soil_ds["bulk_density"] * mask).to_dataframe(name="bulk_density").reset_index().dropna() 

        print(f"Processing {basin} - {crop}")

        # 3 - Baseline Metrics
        sc1_irri = os.path.join(input_dir, "2_2_Sus_Irrigation", f"{basin}_{crop}_annual.csv")
        sc1_rain = os.path.join(input_dir, "2_1_Baseline_rainfed", f"{basin}_{crop}_annual.csv")
        
        # Get metrics + the mass for weighting
        ip_base, ir_base, im_base = get_metrics(sc1_irri, irri_ha, bd_df)
        rp_base, rr_base, rm_base = get_metrics(sc1_rain, rain_ha, bd_df)

        if ip_base is None or rp_base is None:
            continue

        # 4 - Loop Years and Scenarios
        results = []
        for year in range(2009, 2020):
            # Total P Pool = (Pool_irri * Mass_irri + Pool_rain * Mass_rain) / Total_Mass
            total_mass = im_base[year] + rm_base[year]
            base_pool = (ip_base[year] * im_base[year] + rp_base[year] * rm_base[year]) / total_mass
            base_runoff = ir_base[year] + rr_base[year]

            row = {"Basin": basin, "Crop": crop, "Year": year, 
                   "Baseline_P_Pool": base_pool, "Baseline_P_Runoff": base_runoff}

            for sens in Sens_Scenarios:
                i_path = os.path.join(input_dir, "2_3_Sus_Irri_Red_Fert", sens, f"{basin}_{crop}_annual.csv")
                r_path = os.path.join(input_dir, "2_3_Rainfed", sens, f"{basin}_{crop}_annual.csv")
                
                ip_s, ir_s, im_s = get_metrics(i_path, irri_ha, bd_df)
                rp_s, rr_s, rm_s = get_metrics(r_path, rain_ha, bd_df)
                
                if ip_s is not None:
                    row[f"{sens}_P_Pool"] = (ip_s[year] * im_s[year] + rp_s[year] * rm_s[year]) / (im_s[year] + rm_s[year])
                    row[f"{sens}_P_Runoff"] = ir_s[year] + rr_s[year]
            
            results.append(row)

        # 5 - Output
        pd.DataFrame(results).to_csv(os.path.join(output_dir, f"{basin}_{crop}_P_analysis.csv"), index=False)

print("Done.")