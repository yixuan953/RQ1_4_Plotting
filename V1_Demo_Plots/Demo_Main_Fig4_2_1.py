# P runoff over time

import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Studyarea =  ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs/Fig4"

# Define scenarios: "Red_prop" is the main comparison, others are light sensitivity lines
Sens_Scenarios = ["Red_prop", "Red_01", "Red_03", "Red_05", "Red_07", "Red_09", "Red_11", "Red_13", "Red_15"] 

os.makedirs(output_dir, exist_ok=True)

def plot_data(base_p, red_dict, crit_val, crop_name, basin_name):
    plt.figure(figsize=(10, 6))
    
    # 1. Plot Critical Line first (bottom layer)
    plt.axhline(y=1e-6 * crit_val, color='red', linestyle='--', 
                label='Critical P Runoff', linewidth=1.2, zorder=1)
    
    # 2. Plot Sensitivity Lines (Light and Thin)
    # We only add one label for "Sensitivity" to keep the legend clean
    sens_label_added = False
    for folder, red_p in red_dict.items():
        if folder != "Red_prop":
            label = "Sensitivity Scenarios" if not sens_label_added else ""
            plt.plot(red_p.index, 1e-6 * red_p.values, color="#70E8E8", 
                     linewidth=1.2, alpha=0.3, zorder=2, label=label)
            sens_label_added = True

    # 3. Plot Main Reduction Scenario (Proportional)
    if "Red_prop" in red_dict:
        plt.plot(red_dict["Red_prop"].index, 1e-6 * red_dict["Red_prop"].values, 
                 label="Fertilizer Reduction (Main)", color="teal", 
                 marker='o',  linewidth=2.5, zorder=5)
    
    # 4. Plot Baseline (Top layer for visibility)
    plt.plot(base_p.index, 1e-6 * base_p.values, label="Baseline", 
             color="orange", marker='o', linewidth=2.5, zorder=5)
    
    plt.title(f"Annual P Runoff: {crop_name.replace('total ', '').capitalize()} in {basin_name}", fontsize=12)
    plt.xlabel("Year", fontsize=12)
    plt.ylabel("Total P Runoff (ktons P)", fontsize=12)
    
    # Adjust legend: loc='best' or outside if too crowded
    plt.ylim(bottom=0)
    # Ensure all years are visible on x-axis
    plt.xticks(range(2010, 2020), fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='upper left', frameon=True, fontsize=12)
    plt.grid(axis='y', linestyle=':', alpha=0.5)
    plt.tight_layout()
    
    plt.savefig(f"{output_dir}/{basin_name}_{crop_name}_P_runoff.png", dpi=300)
    plt.close()

# --- Main Processing Loop ---

# Re-initialize accumulator for Yangtze rice
yangtze_rice_data = {"baseline": 0, "reduction": {s: 0 for s in Sens_Scenarios}, "critical": 0}

for basin in Studyarea:
    range_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "range.nc")
    with xr.open_dataset(range_nc) as ds_range:
        template = ds_range["mask"]
    low_runoff_path = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)
    
    for crop in Croptypes:
        # Path Logic
        sc1_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline", f"{basin}_{crop}_summary_baseline.nc")
        sc1_irri_path = os.path.join(input_dir, "2_1_Baseline", f"{basin}_{crop}_annual.csv")
        sc1_rain_path = os.path.join(input_dir, "2_1_Baseline_rainfed", f"{basin}_{crop}_annual.csv")

        if not all(os.path.exists(f) for f in [sc1_nc, sc1_irri_path, sc1_rain_path]):
            continue

        # Load Masks
        ds = xr.open_dataset(sc1_nc)
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan) * template * mask_not_low_runoff
        irri_ha_df = (ds["Irrigated_HA"] * mask).to_dataframe(name="Irrigated_HA").reset_index().dropna(subset=["Irrigated_HA"])
        rain_ha_df = (ds["Rainfed_HA"] * mask).to_dataframe(name="Rainfed_HA").reset_index().dropna(subset=["Rainfed_HA"])

        # Calculation Function
        def get_crop_p_runoff(csv_path, area_df, area_col):
            if not os.path.exists(csv_path): return 0
            df = pd.read_csv(csv_path, skipinitialspace=True)
            df = df[(df["Year"] >= 2010) & (df["Year"] <= 2019)]
            merged = df.merge(area_df, left_on=["Lat","Lon"], right_on=["lat","lon"], how="inner")
            merged["P_total"] = merged[area_col] * (merged["P_surf"] + merged["P_sub"])
            return merged.groupby("Year")["P_total"].sum()

        # Calculate Results
        base_p = get_crop_p_runoff(sc1_irri_path, irri_ha_df, "Irrigated_HA") + \
                 get_crop_p_runoff(sc1_rain_path, rain_ha_df, "Rainfed_HA")

        current_crop_red_dict = {}
        for sens in Sens_Scenarios:
            irri_path = os.path.join(input_dir, "2_3_Sus_Irri_Red_Fert", sens, f"{basin}_{crop}_annual.csv")
            rain_path = os.path.join(input_dir, "2_3_Rainfed", sens, f"{basin}_{crop}_annual.csv")
            current_crop_red_dict[sens] = get_crop_p_runoff(irri_path, irri_ha_df, "Irrigated_HA") + \
                                          get_crop_p_runoff(rain_path, rain_ha_df, "Rainfed_HA")

        # Critical Value
        crop_crit_map = {"winterwheat": "Wheat", "maize": "Maize", "soybean": "Soybean", "mainrice": "Rice", "secondrice": "Rice"}
        crit_loss_nc = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/2_Critical_NP_losses/Method3/{crop_crit_map[crop]}/{basin}_crit_P_runoff_kg.nc"
        with xr.open_dataset(crit_loss_nc) as ds_crit:
            crit_val_total = (ds_crit["Critical_P_runoff"] * mask).sum().values.item()

        # Special Case: Yangtze Rice
        if basin == "Yangtze" and (crop == "mainrice" or crop == "secondrice"):
            yangtze_rice_data["baseline"] += base_p
            yangtze_rice_data["critical"] += crit_val_total
            for sens in Sens_Scenarios:
                yangtze_rice_data["reduction"][sens] += current_crop_red_dict[sens]
            
            if crop == "secondrice":
                plot_data(yangtze_rice_data["baseline"], yangtze_rice_data["reduction"], 
                          yangtze_rice_data["critical"], "total rice", "Yangtze")
            continue 

        # Standard Crop Plotting
        plot_data(base_p, current_crop_red_dict, crit_val_total, crop, basin)

