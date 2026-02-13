import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs/Fig4"

os.makedirs(output_dir, exist_ok=True)

def plot_dual_axis(yield_base, yield_red, nue_base, nue_red, crop_name, basin_name):
    plt.rcParams.update({'font.size': 12})
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # --- Left Axis: Yield ---
    color_yield = 'black'
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Total Yield (Mtons)', color=color_yield)
    # Yield plotted as solid lines
    lns1 = ax1.plot(yield_base.index, 1e-9 * yield_base.values, color='orange', marker='o', label='Yield: Baseline', linewidth=2)
    lns2 = ax1.plot(yield_red.index, 1e-9 * yield_red.values, color='teal', marker='o', label='Yield: Fert_Reduction', linewidth=2)
    ax1.tick_params(axis='y', labelcolor=color_yield)
    ax1.set_ylim(bottom=0)

    # --- Right Axis: NUE ---
    ax2 = ax1.twinx()
    color_nue = 'blue'
    ax2.set_ylabel('NUE (Uptake / Total N Input)', color=color_nue)
    # NUE plotted as dashed lines
    lns3 = ax2.plot(nue_base.index, nue_base.values, color='orange', linestyle='--', marker='x', label='NUE: Baseline', alpha=0.7)
    lns4 = ax2.plot(nue_red.index, nue_red.values, color='teal', linestyle='--', marker='x', label='NUE: Fert_Reduction', alpha=0.7)
    ax2.tick_params(axis='y', labelcolor=color_nue)
    ax2.set_ylim(0, 1) # NUE is typically between 0 and 1

    # Combined Legend
    lns = lns1 + lns2 + lns3 + lns4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper left', fontsize=10, frameon=True)

    plt.title(f"Yield & NUE: {crop_name.capitalize()} in {basin_name}")
    plt.xticks(range(2010, 2020))
    ax1.grid(axis='y', linestyle=':', alpha=0.5)
    fig.tight_layout()
    
    plt.savefig(f"{output_dir}/{basin_name}_{crop_name}_Yield_NUE.png", dpi=300)
    plt.close()

# --- Accumulator for Yangtze Rice ---
# Initializing as 0; logic below will handle Series conversion
yangtze_data = {
    "yield_base": 0, "yield_red": 0,
    "uptake_base": 0, "input_base": 0,
    "uptake_red": 0, "input_red": 0
}

for basin in Studyarea:
    range_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "range.nc")
    with xr.open_dataset(range_nc) as ds_range:
        template = ds_range["mask"]
    low_runoff_path = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)

    for crop in Croptypes:
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

        def get_metrics(csv_path, area_df):
            df = pd.read_csv(csv_path, skipinitialspace=True)
            df = df[(df["Year"] >= 2010) & (df["Year"] <= 2019)]
            m = df.merge(area_df, left_on=["Lat","Lon"], right_on=["lat","lon"])
            
            # Totals for Yield and components of NUE
            yield_tot = (m["Storage"] * m["HA"]).groupby(m["Year"]).sum()
            uptake_tot = (m["N_uptake"] * m["HA"]).groupby(m["Year"]).sum()
            
            # Total N Input = Fert + Surf + NH3 + N2O + NOx + N2
            # (Note: Using your specified formula components)
            n_input_tot = (m["HA"] * (m["N_fert"] + m["N_surf"] + m["NH3"] + m["N2O"] + m["NOx"] + m["N2"])).groupby(m["Year"]).sum()
            
            return yield_tot, uptake_tot, n_input_tot

        # Calculate Baseline
        y1_i, u1_i, in1_i = get_metrics(sc1_irri, irri_ha)
        y1_r, u1_r, in1_r = get_metrics(sc1_rain, rain_ha)
        
        y_base = y1_i + y1_r
        nue_base = (u1_i + u1_r) / (in1_i + in1_r)

        # Calculate Reduction
        y2_i, u2_i, in2_i = get_metrics(sc2_irri, irri_ha)
        y2_r, u2_r, in2_r = get_metrics(sc2_rain, rain_ha)
        
        y_red = y2_i + y2_r
        nue_red = (u2_i + u2_r) / (in2_i + in2_r)

        # Yangtze Rice Aggregation
        if basin == "Yangtze" and (crop == "mainrice" or crop == "secondrice"):
            yangtze_data["yield_base"] += y_base
            yangtze_data["yield_red"] += y_red
            yangtze_data["uptake_base"] += (u1_i + u1_r)
            yangtze_data["input_base"] += (in1_i + in1_r)
            yangtze_data["uptake_red"] += (u2_i + u2_r)
            yangtze_data["input_red"] += (in2_i + in2_r)
            
            if crop == "secondrice":
                nue_b_yang = yangtze_data["uptake_base"] / yangtze_data["input_base"]
                nue_r_yang = yangtze_data["uptake_red"] / yangtze_data["input_red"]
                plot_dual_axis(yangtze_data["yield_base"], yangtze_data["yield_red"], 
                               nue_b_yang, nue_r_yang, "total rice", "Yangtze")
            continue

        plot_dual_axis(y_base, y_red, nue_base, nue_red, crop, basin)