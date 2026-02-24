# This script is used to plot the P pool changes over time in different scenarios

import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configuration
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
CropTypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
# Define scenarios: "Red_prop" is the main comparison, others are light sensitivity lines
Sens_Scenarios = ["Red_prop", "Red_01", "Red_03", "Red_05", "Red_07", "Red_09", "Red_11", "Red_13", "Red_15"] 

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs/Fig4"
os.makedirs(output_dir, exist_ok=True)

p_mol = 30.97 # 1 mol P = 30.97 g
crit_P_olsen = 15 # [mg/kg]
Trans_fact = 8 # P Pxalate = P Olsen * 8 = Labile + Stable P pool

def get_metrics(csv_path, area_df, bd_df):
    df = pd.read_csv(csv_path, skipinitialspace=True)
    # Standardize column names to lowercase to match the .nc dataframes
    df = df.rename(columns={'Lat': 'lat', 'Lon': 'lon'})
    
    df = df[(df["Year"] >= 2010) & (df["Year"] <= 2019)]
    
    # Merge datasets on lowercase lat/lon
    m = df.merge(area_df, on=["lat", "lon"]).merge(bd_df, on=["lat", "lon"])
    
    # Calculate weighted components
    m["pixel_mass"] = m["HA"] * m["bulk_density"]
    m["weighted_P"] = (m["LabileP"] + m["StableP"]) * m["pixel_mass"] * p_mol
    
    # Group by Year
    grouped = m.groupby("Year")

    p_pool = grouped["weighted_P"].sum() / grouped["pixel_mass"].sum()
    
    return p_pool

def plot_data(base_p, red_dict, crop_name, basin_name, irrigation):
    plt.figure(figsize=(10, 6))

    # 1. Plot the Grey Range (shading)
    # We calculate the lower and upper bounds using the +/- 3 margin
    lower_bound = (crit_P_olsen - 3) * Trans_fact
    upper_bound = (crit_P_olsen + 3) * Trans_fact

    plt.axhspan(lower_bound, upper_bound, color='grey', alpha=0.3, 
                label='P pool size for crop growth without P limitation\n[P Olsen = 15 ± 3 mg/kg]', zorder=0)
    plt.axhline(y=Trans_fact * crit_P_olsen, color='red', linestyle='--', 
                label='Critical N Runoff', linewidth=1.2, zorder=1)

    # 2. Plot Sensitivity Lines (Light and Thin)
    # We only add one label for "Sensitivity" to keep the legend clean
    sens_label_added = False
    for folder, red_p in red_dict.items():
        if folder != "Red_prop":
            label = "Sensitivity Scenarios" if not sens_label_added else ""
            plt.plot(red_p.index, red_p.values, color="#70E8E8", 
                     linewidth=1.2, alpha=0.3, zorder=2, label=label)
            sens_label_added = True

    # 3. Plot Main Reduction Scenario (Proportional)
    if "Red_prop" in red_dict:
        plt.plot(red_dict["Red_prop"].index, red_dict["Red_prop"].values, 
                 label="Fertilizer Reduction (Main)", color="teal", 
                 marker='o',  linewidth=2.5, zorder=5)
    
    # 4. Plot Baseline (Top layer for visibility)
    plt.plot(base_p.index,  base_p.values, label="Baseline", 
             color="orange", marker='o', linewidth=2.5, zorder=5)
    
    plt.title(f"Topsoil (30 cm) P pool changes: {crop_name.replace('total ', '').capitalize()} in {basin_name}", fontsize=12)
    plt.xlabel("Year", fontsize=12)
    plt.ylabel("Labile + Stable P pool [mg/kg]", fontsize=12)
    
    # Adjust legend: loc='best' or outside if too crowded
    plt.ylim(0, 400)
    # Ensure all years are visible on x-axis
    plt.xticks(range(2010, 2020), fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='upper left', frameon=True, fontsize=12)
    plt.grid(axis='y', linestyle=':', alpha=0.5)
    plt.tight_layout()
    
    plt.savefig(f"{output_dir}/{basin_name}_{crop_name}_soil_P_pool_{irrigation}.png", dpi=300)
    plt.close()

for basin in Studyareas:
    range_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "range.nc")


    with xr.open_dataset(range_nc) as ds_range:
        template = ds_range["mask"]

    low_runoff_path = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan) 

    for crop in CropTypes:

        # 1 - Paths for Harvested area
        sc1_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline", f"{basin}_{crop}_summary_baseline.nc")
        if not os.path.exists(sc1_nc):
            print (f"Skipping {basin} - {crop}")
            continue
        ds = xr.open_dataset(sc1_nc)
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan) * template * mask_not_low_runoff
        irri_ha = (ds["Irrigated_HA"] * mask).to_dataframe(name="HA").reset_index().dropna()
        rain_ha = (ds["Rainfed_HA"] * mask).to_dataframe(name="HA").reset_index().dropna()

        # 2 - Path for soil bulk density file
        soil_prop_file = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "Mask", f"{basin}_{crop}_mask.nc")
        if not os.path.exists(soil_prop_file):
            continue
        soil_ds = xr.open_dataset(soil_prop_file)
        bd_df = (soil_ds["bulk_density"] * mask).to_dataframe(name="bulk_density").reset_index().dropna() 

        print(f"Processing {basin} - {crop}")

        # 3 - Paths model outputs
        # Scenario 1 - Baseline
        sc1_irri = os.path.join(input_dir, "2_1_Baseline", f"{basin}_{crop}_annual.csv")
        sc1_rain = os.path.join(input_dir, "2_1_Baseline_rainfed", f"{basin}_{crop}_annual.csv")
        if not os.path.exists(sc1_irri) or not os.path.exists(sc1_rain):
            print (f"{basin} does not have {crop} baseline .csv")
            continue
        Irri_base_P_pool = get_metrics(sc1_irri, irri_ha, bd_df)
        Rainfed_base_P_pool = get_metrics(sc1_rain, rain_ha, bd_df)

        # Scenario 2 - Reduced fertilization
        Irrigated_crop_red_dict = {}
        Rainfed_crop_red_dict = {}

        for sens in Sens_Scenarios:
            irri_path = os.path.join(input_dir, "2_3_Sus_Irri_Red_Fert", sens, f"{basin}_{crop}_annual.csv")
            rain_path = os.path.join(input_dir, "2_3_Rainfed", sens, f"{basin}_{crop}_annual.csv")

            Irrigated_crop_red_dict[sens] = get_metrics(irri_path, irri_ha, bd_df)
            Rainfed_crop_red_dict[sens] = get_metrics(rain_path, rain_ha, bd_df)
            
        
        plot_data(Irri_base_P_pool, Irrigated_crop_red_dict, crop, basin, irrigation = "Irri")
        plot_data(Rainfed_base_P_pool, Rainfed_crop_red_dict, crop, basin, irrigation = "Rain")

