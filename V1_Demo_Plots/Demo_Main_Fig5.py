# This script is used to prepare the annual: Fertilizer input, NUE, P pool size, and transpiration

import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configuration
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
CropTypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs/Fig5"
os.makedirs(output_dir, exist_ok=True)

def get_metrics(csv_path, area_df, bd_df):
    df = pd.read_csv(csv_path, skipinitialspace=True)
    df = df[(df["Year"] >= 2010) & (df["Year"] <= 2019)]
    m1 = df.merge(area_df, left_on=["Lat","Lon"], right_on=["lat","lon"])
    m = m1.merge(bd_df, left_on=["Lat","Lon"], right_on=["lat","lon"])
            
    # Interested variables
    n_input_tot = (m["HA"] * (m["N_fert"] + m["N_surf"] + m["NH3"] + m["N2O"] + m["NOx"] + m["N2"])).groupby(m["Year"]).sum()
    n_uptake_tot = (m["N_uptake"] * m["HA"]).groupby(m["Year"]).sum() # kg
    p_uptake_tot = (m["P_uptake"] * m["HA"]).groupby(m["Year"]).sum() # kg
    p_acc = (m["P_acc"] * m["HA"]).groupby(m["Year"]).sum()           # kg
    pool_den = (m["HA"] * m["bulk_density"]).sum()
    p_pool = ((m["LabileP"] + m["StableP"]) * m["HA"] * m["bulk_density"]).groupby(m["Year"]).sum()/pool_den # mmol/kg
    np_uptake_ratio = n_uptake_tot/p_uptake_tot 
            
    return n_input_tot, n_uptake_tot, p_uptake_tot, p_acc, p_pool, np_uptake_ratio

def plot(n_input, n_uptake, p_uptake, p_acc, p_pool, np_ratio, 
         n_input_red, n_uptake_red, p_uptake_red, p_acc_red, p_pool_red, np_ratio_red, 
         crop_name, basin_name, plot_type):

    # 1. Prepare Data: Convert units and calculate NUE series
    # (kg to ktons: / 1e6)
    data_baseline = [
        n_input / 1e6, 
        n_uptake / 1e6, 
        n_uptake / n_input, # NUE
        p_pool, 
        p_uptake / 1e6, 
        p_acc / 1e6,
        np_ratio
    ]
    
    data_red = [
        n_input_red / 1e6, 
        n_uptake_red / 1e6, 
        n_uptake_red / n_input_red, # NUE
        p_pool_red, 
        p_uptake_red / 1e6, 
        p_acc_red / 1e6,
        np_ratio_red
    ]

    labels = ['N Input\n(ktons)', 'N Uptake\n(ktons)', 'NUE', 
              'P Pool\n(mmol/kg)', 'P Uptake\n(ktons)', 'P Accum.\n(ktons)', 'N:P uptake']

    # 2. Setup Figure
    fig, axes = plt.subplots(7, 1, figsize=(3.5, 16), sharex=True)
    fig.subplots_adjust(hspace=0.5)
    
    # Aesthetic colors
    c_base = "#f4ae64"  # Soft Salmon
    c_red = "#4bbd8f" # Sage Green

    for i in range(7):
        ax = axes[i]
        plot_data = [data_baseline[i], data_red[i]]
        positions = [1.4, 1.8]
        colors = [c_base, c_red]

        for pos, d, color in zip(positions, plot_data, colors):
            # 1. Calculate stats
            median = d.median()
            q1 = d.quantile(0.25)
            q3 = d.quantile(0.75)
                
            # 2. Draw the "Stem" (Interquartile Range) with rounded caps
            ax.plot([pos, pos], [q1, q3], color=color, linewidth=6, solid_capstyle='round', alpha=0.4)
                
            # 3. Draw the "Whisker" (Full range) - very thin
            ax.plot([pos, pos], [d.min(), d.max()], color=color, linewidth=1, linestyle='-', alpha=0.6)

            # 4. Draw the "Bubble" (Median) - The cute circular part
            ax.scatter(pos, median, s=120, color=color, edgecolor='white', linewidth=2, zorder=5)

            # 5. The "Seeds" (Jittered Year Points)
            x_jitter = np.random.normal(pos, 0.05, size=len(d))
            ax.scatter(x_jitter, d, s=15, color='#2c3e50', alpha=0.3, edgecolor='none', zorder=4)

        # Styling
        ax.set_ylabel(labels[i], fontsize=9, fontweight='bold', color='#4f4f4f')
        ax.set_xticks([1.4, 1.8])
        ax.set_xticklabels(['Baseline', 'Fert_Red'], fontsize=10, fontweight='bold')
        
        # Clean up axes
        for spine in ['top', 'right', 'bottom']:
            ax.spines[spine].set_visible(False)
        ax.spines['left'].set_color('#eeeeee')
        ax.tick_params(axis='both', which='both', length=0) # Hide tick marks for a cleaner look
        ax.grid(axis='y', linestyle='-', color='#f0f0f0', zorder=0)

        # Styling
        ax.set_ylabel(labels[i], fontsize=9, fontweight='bold', color='#4f4f4f')
        ax.set_xticks([1.4, 1.8])
        ax.set_xticklabels(['Baseline', 'Reduction'], fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.3)

    plt.suptitle(f"{basin_name.upper()} - {crop_name.capitalize()} ({plot_type})\nAnnual Variability (2010-2019)", 
                 fontsize=14, fontweight='black', y=0.93)

    # 3. Save with plot_type in name to avoid overwriting
    output_path = f"{output_dir}/{basin_name}_{crop_name}_{plot_type}_boxplot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
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
        soil_prop_file = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "Mask", f"{basin}_{crop}_mask.nc")
        if not os.path.exists(soil_prop_file):
            continue
        print(f"Processing {basin} - {crop}")
    
        # Paths for Scenario 1 (Baseline) and Scenario 2 (Main Reduction - Red_prop)
        sc1_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline", f"{basin}_{crop}_summary_baseline.nc")
        
        # Scenario 1 Files
        sc1_irri = os.path.join(input_dir, "2_1_Baseline", f"{basin}_{crop}_annual.csv")
        sc1_rain = os.path.join(input_dir, "2_1_Baseline_rainfed", f"{basin}_{crop}_annual.csv")
        
        # Scenario 2 Files (Main Reduction)
        sc2_irri = os.path.join(input_dir, "2_3_Sus_Irri_Red_Fert", "Red_prop", f"{basin}_{crop}_annual.csv")
        sc2_rain = os.path.join(input_dir, "2_3_Rainfed", "Red_prop", f"{basin}_{crop}_annual.csv")

        if not all(os.path.exists(f) for f in [sc1_nc, sc1_irri, sc1_rain, sc2_irri, sc2_rain]):
            continue

        ds = xr.open_dataset(sc1_nc)
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan) * template * mask_not_low_runoff
        irri_ha = (ds["Irrigated_HA"] * mask).to_dataframe(name="HA").reset_index().dropna()
        rain_ha = (ds["Rainfed_HA"] * mask).to_dataframe(name="HA").reset_index().dropna()

        soil_ds = xr.open_dataset(soil_prop_file)
        bd_df = (soil_ds["bulk_density"] * mask).to_dataframe(name="bulk_density").reset_index().dropna()

        Irri_n_input_tot_baseline, Irri_n_uptake_tot_baseline, Irri_p_uptake_tot_baseline, Irri_p_acc_baseline, Irri_p_pool_baseline, Irri_np_uptake_ratio_baseline  = get_metrics(sc1_irri, irri_ha, bd_df)
        Irri_n_input_tot_red, Irri_n_uptake_tot_red, Irri_p_uptake_tot_red, Irri_p_acc_red, Irri_p_pool_red, Irri_np_uptake_ratio_red = get_metrics(sc2_irri, irri_ha, bd_df)

        Rain_n_input_tot_baseline, Rain_n_uptake_tot_baseline, Rain_p_uptake_tot_baseline, Rain_p_acc_baseline, Rain_p_pool_baseline, Rain_np_uptake_ratio_baseline = get_metrics(sc1_rain, rain_ha, bd_df)
        Rain_n_input_tot_red, Rain_n_uptake_tot_red, Rain_p_uptake_tot_red, Rain_p_acc_red, Rain_p_pool_red, Rain_np_uptake_ratio_red = get_metrics(sc2_rain, rain_ha, bd_df)

        plot(Irri_n_input_tot_baseline, Irri_n_uptake_tot_baseline, Irri_p_uptake_tot_baseline, Irri_p_acc_baseline, Irri_p_pool_baseline, Irri_np_uptake_ratio_baseline, Irri_n_input_tot_red, Irri_n_uptake_tot_red, Irri_p_uptake_tot_red, Irri_p_acc_red, Irri_p_pool_red, Irri_np_uptake_ratio_red, crop, basin, plot_type="Irrigated")
        plot(Rain_n_input_tot_baseline, Rain_n_uptake_tot_baseline, Rain_p_uptake_tot_baseline, Rain_p_acc_baseline, Rain_p_pool_baseline, Rain_np_uptake_ratio_baseline, Rain_n_input_tot_red, Rain_n_uptake_tot_red, Rain_p_uptake_tot_red, Rain_p_acc_red, Rain_p_pool_red, Rain_np_uptake_ratio_red, crop, basin, plot_type="Rainfed")