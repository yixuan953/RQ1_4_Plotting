# This script is used to plot the 2D exceedance of N and P boundary exceedance

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import cartopy.crs as ccrs
import geopandas as gpd
import pandas as pd


DATA_DIR = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
BASE_DIR = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
SCEN_DIR = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
FIG_OUT  = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig4_Crop_Prod_Red"
STUDY_AREAS = ["LaPlata", "Rhine", "Indus", "Yangtze"]
CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}


def get_bivariate_data(folder, basin, crop_list):
    """Aggregates multiple sub-crops and calculates exceedance."""
    n_run_total = None
    p_run_total = None
    n_crit_total = None
    p_crit_total = None
    
    lr_path = os.path.join(DATA_DIR, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(lr_path) as ds_lr:
        low_flow_mask = ds_lr["Low_Runoff"]
        non_low_flow_area = xr.where(ds_lr["Low_Runoff"] == 1, np.nan, 1)

    data_found = False
    for crop in crop_list:
        path = os.path.join(folder, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(path): 
            continue
        
        data_found = True
        with xr.open_dataset(path) as ds:
            basin_mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
            mask = basin_mask * non_low_flow_area
            
            n_r = ds["N_Runoff"] * mask
            p_r = ds["P_Runoff"] * mask
            n_c = ds["Crit_N_Runoff"] * mask
            p_c = ds["Crit_P_Runoff"] * mask

            if n_run_total is None:
                n_run_total, p_run_total = n_r, p_r
                n_crit_total, p_crit_total = n_c, p_c
            else:
                n_run_total = n_run_total.fillna(0) + n_r.fillna(0)
                p_run_total = p_run_total.fillna(0) + p_r.fillna(0)
                n_crit_total = n_crit_total.fillna(0) + n_c.fillna(0)
                p_crit_total = p_crit_total.fillna(0) + p_c.fillna(0)

    # If no files were found for this basin/crop combination, return None
    if not data_found:
        return None, None, low_flow_mask

    # Calculate final exceedance
    n_exc = 100 * (n_run_total - n_crit_total) / n_crit_total
    p_exc = 100 * (p_run_total - p_crit_total) / p_crit_total
    
    return n_exc, p_exc, low_flow_mask

def plot_bivariate_grid():
    # 4x4 Bivariate Palette
    colors = [
        # N-low (Low N/P is now a clearer Mint Green)
        "#e2f0d9", "#f0cc89", "#f1b346", "#ef9b00", # Pure P (Row 1) remains Gold/Orange
        
        # N-inc (Brighter Orchid)
        "#f48fb1", "#d4927d", "#c5794f", "#b65f21", 
        
        # N-high (Vibrant Rose)
        "#f06292", "#b26842", "#a04020", "#9e4216", 
        
        # N-very high (The "Electric Rose" #E41F7E -> Deep Brown)
        "#E41F7E", "#8c4a2f", "#7a301a", "#512611"  # Pure N (Bottom Left) is #E41F7E
    ]
    cmap_2d = mcolors.ListedColormap(colors)

    crops = list(CROP_MAP.keys())
    fig, axes = plt.subplots(len(crops), len(STUDY_AREAS), figsize=(24, 20), 
                             subplot_kw={'projection': ccrs.PlateCarree()})

    for c, basin in enumerate(STUDY_AREAS):
        shp_path = os.path.join(DATA_DIR, "2_shp_StudyArea", basin, f"{basin}.shp")
        shp = gpd.read_file(shp_path)
        bbox = shp.total_bounds

        for r, crop_name in enumerate(crops):
            ax = axes[r, c]
            ax.axis('off')
            
            n_exc, p_exc, low_flow = get_bivariate_data(SCEN_DIR, basin, CROP_MAP[crop_name])
            
            if n_exc is not None:
                # 1. Bivariate Index Calculation
                def bin_data(da):
                    # Clips values < 0 to 0 so they fall into the first bin
                    # Any exceedance <= 25 (including negative) becomes bin 0
                    return xr.where(da <= 25, 0, 
                           xr.where(da <= 50, 1, 
                           xr.where(da <= 75, 2, 3)))

                # Calculate the 0-15 index
                # If both are <= 25% (including < 0%), the result is 0 (Mint Green)
                idx = (bin_data(n_exc) * 4 + bin_data(p_exc))
                
                # IMPORTANT: Only mask the 'No Data' areas (where n_exc itself is NaN)
                # This keeps the 'safe' (exceedance < 0) pixels as index 0
                idx = idx.where(n_exc.notnull())

                idx.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap_2d, 
                         add_colorbar=False, vmin=0, vmax=15, zorder=1)

                # 2. Masked Grey Low Flow
                # Ensure the grey only covers pixels that have valid crop data
                masked_low_flow = low_flow
                
                masked_low_flow.plot(ax=ax, transform=ccrs.PlateCarree(), 
                                     cmap=mcolors.ListedColormap(["#e7e7e7"]), 
                                     add_colorbar=False, zorder=2)

                # 3. Basin Boundary - Top Layer
                shp.boundary.plot(ax=ax, color='black', linewidth=0.8, zorder=3)
                ax.set_extent([bbox[0], bbox[2], bbox[1], bbox[3]], crs=ccrs.PlateCarree())
            else:
                ax.text(0.5, 0.5, "No Data", transform=ax.transAxes, 
                        ha='center', va='center', color='grey', fontsize=18)

    # Detail 1: 2D Bivariate Legend with 0-100 range
    ax_leg = fig.add_axes([0.42, 0.04, 0.12, 0.12]) 
    leg_data = np.arange(16).reshape(4, 4)
    ax_leg.imshow(leg_data, cmap=cmap_2d, origin='lower', extent=[0, 100, 0, 100])
    
    ax_leg.set_xlabel("P [%]", fontsize=20)
    ax_leg.set_ylabel("N [%]", fontsize=20)
    
    # Set ticks to show the 0-100 scale
    ticks = [0, 25, 50, 75, 100]
    ax_leg.set_xticks(ticks)
    ax_leg.set_yticks(ticks)
    ax_leg.tick_params(labelsize=18)

    plt.savefig(os.path.join(FIG_OUT, "Bivariate_NP_Final_Grid.png"), dpi=300, bbox_inches='tight')

plot_bivariate_grid()


def analyze_extreme_exceedance(threshold=50):
    """
    Calculates the proportion of Harvested Area (HA) where 
    N or P exceedance is > threshold.
    """
    stats_records = []

    for basin in STUDY_AREAS:
        for crop_name, crop_list in CROP_MAP.items():
            # Use your existing data fetching logic
            n_exc, p_exc, _ = get_bivariate_data(SCEN_DIR, basin, crop_list)
            
            if n_exc is None:
                continue

            # Load the Total Harvested Area for these crops
            # We sum HA across the crop_list (e.g., mainrice + secondrice)
            total_ha_basin_crop = 0
            extreme_n_ha = 0
            extreme_p_ha = 0
            extreme_both_ha = 0

            for crop in crop_list:
                path = os.path.join(SCEN_DIR, f"{basin}_{crop}_summary.nc")
                if not os.path.exists(path): continue
                
                with xr.open_dataset(path) as ds:
                    ha = ds["Total_HA"].where(ds["Total_HA"] > 2500, 0)
                    
                    # Calculate exceedance for this specific sub-crop
                    n_e = 100 * (ds["N_Runoff"] - ds["Crit_N_Runoff"]) / ds["Crit_N_Runoff"]
                    p_e = 100 * (ds["P_Runoff"] - ds["Crit_P_Runoff"]) / ds["Crit_P_Runoff"]

                    # Sum up the Area (HA) meeting the "extreme" criteria
                    total_ha_basin_crop += ha.sum().values
                    extreme_n_ha += ha.where(n_e > threshold).sum().values
                    extreme_p_ha += ha.where(p_e > threshold).sum().values
                    extreme_both_ha += ha.where((n_e > threshold) & (p_e > threshold)).sum().values

            # Calculate Proportions (%)
            if total_ha_basin_crop > 0:
                stats_records.append({
                    "Basin": basin,
                    "Crop": crop_name,
                    "Total_HA_kHa": total_ha_basin_crop / 1000,
                    "Extreme_N_%": (extreme_n_ha / total_ha_basin_crop) * 100,
                    "Extreme_P_%": (extreme_p_ha / total_ha_basin_crop) * 100,
                    "Extreme_Both_%": (extreme_both_ha / total_ha_basin_crop) * 100
                })

    df_stats = pd.DataFrame(stats_records)
    return df_stats

# Run the analysis
df_results = analyze_extreme_exceedance(threshold=100)
print(df_results.to_string(index=False))

# Optional: Save to CSV for your report
df_results.to_csv(os.path.join(FIG_OUT, "Extreme_Exceedance_Stats.csv"))