# This code is used to plot: 
# 1 - Rainfed exceedance [kg N/ha]
# 2 - Irrigated exceedance [kg N/ha]
# 3 - Total exceedance [kg N/ha]

# 4 - Rainfed exceedance [kg P/ha]
# 5 - Irrigated exceedance [kg P/ha]
# 6 - Total exceedance [kg N/ha]

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap
import cartopy.geodesic as cgeo

# --- Configuration ---
Studyareas = ["Yangtze"] # ["LaPlata", "Indus", "Yangtze", "Rhine"]
CropGroups = ["winterwheat", "maize", "soybean", "Rice"] 

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig3_Exceedance"

grey_cmap = ListedColormap(["#EBE9E9"])

# Uniform ranges for kg/ha across N and P
VAR_RANGES = {
    "N_Runoff_kg_ha": (-20, 20),
    "P_Runoff_kg_ha": (-1.0, 1.0),
}

def GetSimSummary(ds):
    """Calculates exceedance variables in kg/ha only."""
    mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)

    # 1. BOUNDARY RATE [kg/ha]
    # Critical Load [kg] / Total Area [ha]
    Crit_N_rate = (ds["Crit_N_Runoff"] / ds["Total_HA"]).where(ds["Total_HA"] > 0) * mask
    Crit_P_rate = (ds["Crit_P_Runoff"] / ds["Total_HA"]).where(ds["Total_HA"] > 0) * mask
    
    # 2. CURRENT RATES [kg/ha]
    # Rainfed and Irrigated are already in kg/ha in the .nc
    N_rain = ds["N_Runoff_Rainfed"].where(ds["Rainfed_HA"] > 250) * mask
    N_irri = ds["N_Runoff_Irrigated"].where(ds["Irrigated_HA"] > 250) * mask
    N_avg  = ds["N_Runoff"]/ds["Total_HA"].where(ds["Total_HA"] > 0)  * mask # This is the basin-wide average rate

    P_rain = ds["P_Runoff_Rainfed"].where(ds["Rainfed_HA"] > 250) * mask
    P_irri = ds["P_Runoff_Irrigated"].where(ds["Irrigated_HA"] > 250) * mask
    P_avg  = ds["P_Runoff"]/ds["Total_HA"].where(ds["Total_HA"] > 0)  * mask

    # 3. EXCEEDANCE [Current - Boundary]
    return (mask, 
            N_rain - Crit_N_rate,  # 1: Rainfed N Exceed [kg/ha]
            N_irri - Crit_N_rate,  # 2: Irrigated N Exceed [kg/ha]
            N_avg  - Crit_N_rate,  # 3: Average N Exceed [kg/ha]
            P_rain - Crit_P_rate,  # 4: Rainfed P Exceed [kg/ha]
            P_irri - Crit_P_rate,  # 5: Irrigated P Exceed [kg/ha]
            P_avg  - Crit_P_rate   # 6: Average P Exceed [kg/ha]
           )

def calculate_stats(data_array):
    vals = data_array.values.flatten()
    vals = vals[~np.isnan(vals) & (vals != 0)]
    if len(vals) == 0: return np.nan, np.nan, np.nan
    return np.percentile(vals, 25), np.percentile(vals, 50), np.percentile(vals, 75)

def add_scale_bar(ax, length_km, basin, location_y=-0.08):
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    x_start = lon0 + (lon1 - lon0) * 0.05
    y_pos = lat0 + (lat1 - lat0) * location_y
    ax.plot([x_start, x_start + bar_width_deg], [y_pos, y_pos], transform=ccrs.PlateCarree(),
            color='black', linewidth=8, zorder=10, clip_on=False)
    ax.text(x_start, y_pos - (lat1 - lat0) * 0.03, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha='left', va='top', fontsize=70, clip_on=False)

def plot_boundaries(ds, data_dir, fig_path, basin):
    outputs = GetSimSummary(ds)
    mask = outputs[0]
    data_list = outputs[1:] 
    
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        lr_area = ds_lr["Low_Runoff"].where(mask.notnull())

    # Map subplots to range keys
    range_map = ["N_Runoff_kg_ha"]*3 + ["P_Runoff_kg_ha"]*3
    titles = ["Rainfed N", "Irrigated N", "Average N", "Rainfed P", "Irrigated P", "Average P"]
    
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    gdf_boundary = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    basin_aspect = (lat_max - lat_min) / (lon_max - lon_min)

    fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(30, 130),
                             subplot_kw={'projection': ccrs.PlateCarree()})
    
    for i, ax in enumerate(axes):
        vmin, vmax = VAR_RANGES[range_map[i]]
        ax.set_box_aspect(basin_aspect)
        
        im = data_list[i].plot(ax=ax, transform=ccrs.PlateCarree(), 
                               cmap='RdYlGn_r', vmin=vmin, vmax=vmax, add_colorbar=False)
        
        lr_area.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=grey_cmap, add_colorbar=False, zorder=2)
        gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=4.0, zorder=3)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.axis('off')
        # ax.text(0.02, 0.95, titles[i], transform=ax.transAxes, fontsize=80, weight='bold', va='top')

        # # --- Colorbar ---
        # pos = ax.get_position()
        # cax = fig.add_axes([pos.x0, pos.y0 - 0.04, pos.width, 0.012])
        # cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
        # cbar.ax.tick_params(labelsize=80, length=20, width=5)

        if i == 0:  
            bar_len = 500 if (lon_max - lon_min) > 10 else 100
            add_scale_bar(ax, bar_len, basin)

    plt.subplots_adjust(hspace=0.85) 
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

# --- Loop ---
stats_results = []
for basin in Studyareas:
    for group in CropGroups:
        sub_crops = ["secondrice","mainrice"] if (group == "Rice" and basin == "Yangtze") else \
                    (["mainrice"] if group == "Rice" else [group])

        datasets = [xr.open_dataset(os.path.join(input_dir, f"{basin}_{sc}_summary.nc")) 
                    for sc in sub_crops if os.path.exists(os.path.join(input_dir, f"{basin}_{sc}_summary.nc"))]
        
        if not datasets: continue
        
        if len(datasets) > 1:
            # 1. Use the first dataset as a base, but fill its NaNs with 0 
            # to ensure we don't 'kill' data from the second dataset during the sum
            ds_combined = datasets[0].fillna(0).copy()
            
            # 2. Sum the Area variables [ha]
            # .fillna(0) is critical here: it ensures that (MainRice + 0) = MainRice
            area_vars = ["Total_HA", "Rainfed_HA", "Irrigated_HA"]
            for v in area_vars:
                ds_combined[v] = sum(d[v].fillna(0) for d in datasets)
            
            # 3. Sum the Critical Load Mass variables [kg]
            mass_vars = ["Crit_N_Runoff", "Crit_P_Runoff"]
            for v in mass_vars:
                ds_combined[v] = sum(d[v].fillna(0) for d in datasets)

            # 4. Area-Weighted Average for Intensities [kg/ha]
            # We calculate Total Mass / Total Area for each pixel
            for v_base in ["N_Runoff", "P_Runoff"]:
                # Total Average
                total_mass = sum((d[v_base] * d["Total_HA"]).fillna(0) for d in datasets)
                ds_combined[v_base] = total_mass / ds_combined["Total_HA"].where(ds_combined["Total_HA"] > 0)
                
                # Rainfed Average
                rain_mass = sum((d[f"{v_base}_Rainfed"] * d["Rainfed_HA"]).fillna(0) for d in datasets)
                ds_combined[f"{v_base}_Rainfed"] = rain_mass / ds_combined["Rainfed_HA"].where(ds_combined["Rainfed_HA"] > 0)
                
                # Irrigated Average
                irri_mass = sum((d[f"{v_base}_Irrigated"] * d["Irrigated_HA"]).fillna(0) for d in datasets)
                ds_combined[f"{v_base}_Irrigated"] = irri_mass / ds_combined["Irrigated_HA"].where(ds_combined["Irrigated_HA"] > 0)

            # 5. Fix the Basin Mask
            # Ensure the mask is 1 wherever ANY of the crops exist
            ds_combined["Basin_mask"] = sum(d["Basin_mask"].fillna(0) for d in datasets).clip(max=1)
            
        else:
            ds_combined = datasets[0]

        output_file = os.path.join(fig_base_dir, f"{basin}_{group}_exceedance_kgha.png")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plot_boundaries(ds_combined, data_dir, output_file, basin)

        # Statistics
        outputs = GetSimSummary(ds_combined)
        vars = ["N_Rain_ex", "N_Irri_ex", "N_Avg_ex", "P_Rain_ex", "P_Irri_ex", "P_Avg_ex"]
        for idx, v_name in enumerate(vars):
            q1, med, q3 = calculate_stats(outputs[idx+1])
            stats_results.append({"Basin": basin, "Crop": group, "Variable": v_name, "25th": q1, "Median": med, "75th": q3})

pd.DataFrame(stats_results).to_csv(os.path.join(fig_base_dir, "Summary_Stats_kgha.csv"), index=False)
print("Done.")