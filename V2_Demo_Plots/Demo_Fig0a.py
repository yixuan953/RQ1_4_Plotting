# This code is used to plot: 
# 1 - Rainfed yield [kg/ha]
# 2 - Irrigated yield [kg/ha]
# 3 - Total crop production [kg/ha]

# 4 - Rainfed N runoff [kg/ha
# 5 - Irrigated N runoff [kg/ha]
# 6 - Critical N runoff [kg/ha]
# 7 - Total N runoff [ktons]
# 8 - Critical N runoff [ktons]

# 9 - Rainfed P runoff [kg/ha]
# 10 - Irrigated P runoff [kg/ha]
# 11 - Critical P runoff [kg/ha]
# 12 - Total P runoff [ktons]
# 13 - Critical P runoff [ktons]

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
from matplotlib.colors import LinearSegmentedColormap
import cartopy.geodesic as cgeo
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# --- Configuration ---
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
InputCrops = ["winterwheat"]
# InputCrops = ["mainrice", "secondrice"]  # Uncomment to include all crops
# InputCrops = ["soybean"]  # Uncomment to include all crops
# InputCrops = ["maize"]  # Uncomment to include all crops

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig0"

# Custom Phosphorus Colormap
p_colors = ["#fff6e9f1", "#feb728"]
custom_phosphorus_cmap = LinearSegmentedColormap.from_list("custom_p", p_colors, N=256)

grey_cmap = ListedColormap(["#EBE9E9"]) # Standard Grey


VAR_RANGES = {
    "Yield_kg_ha":      (0, 6000),
    "Production_ktons": (0, 200),  
    "N_Runoff_kg_ha":   (0, 80),
    "N_Runoff_ktons":   (0, 2),
    "P_Runoff_kg_ha":   (0, 3),
    "P_Runoff_ktons":   (0, 0.1)
}

def GetSimSummary(file_name):
    ds = xr.open_dataset(file_name)
    # Threshold masking
    mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)

    Irrigated_Yield = ds["Avg_Yield_Irrigated"].where(ds["Irrigated_HA"] > 250) * mask
    Rainfed_Yield = ds["Avg_Yield_Rainfed"].where(ds["Rainfed_HA"] > 250) * mask
    
    # Production calculation (kg to ktons)
    Total_Production = 1e-6 * mask * (
        ds["Avg_Yield_Irrigated"].fillna(0) * ds["Irrigated_HA"].fillna(0) + 
        ds["Avg_Yield_Rainfed"].fillna(0) * ds["Rainfed_HA"].fillna(0)
    )

    # N Runoff logic (Note: Fix indices/names if they differ in your .nc)
    N_runoff_Rainfed   = ds["N_Runoff_Rainfed"].where(ds["Rainfed_HA"] > 250) * mask
    N_runoff_Irrigated = ds["N_Runoff_Irrigated"].where(ds["Irrigated_HA"] > 250) * mask
    Crit_N_kgperha     = (ds["Crit_N_Runoff"] / ds["Total_HA"]) * mask

    N_runoff = 1e-6 * ds["N_Runoff"] * mask
    Crit_N   = 1e-6 * ds["Crit_N_Runoff"] * mask

    # P Runoff logic
    P_runoff_Rainfed   = ds["P_Runoff_Rainfed"].where(ds["Rainfed_HA"] > 250) * mask
    P_runoff_Irrigated = ds["P_Runoff_Irrigated"].where(ds["Irrigated_HA"] > 250) * mask
    Crit_P_kgperha     = (ds["Crit_P_Runoff"] / ds["Total_HA"]) * mask

    P_runoff = 1e-6 * ds["P_Runoff"] * mask
    Crit_P   = 1e-6 * ds["Crit_P_Runoff"] * mask
    
    return (mask, Rainfed_Yield, Irrigated_Yield, Total_Production, 
            N_runoff_Rainfed, N_runoff_Irrigated, Crit_N_kgperha, N_runoff, Crit_N, 
            P_runoff_Rainfed, P_runoff_Irrigated, Crit_P_kgperha, P_runoff, Crit_P)


def add_scale_bar(ax, length_km, basin, location_y=-0.15):
    """
    Places the scale bar below the plot.
    """
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    
    if basin in ["Yangtze", "Rhine"]:
        x_start = lon0 + (lon1 - lon0) * 0.05
        x_end = x_start + bar_width_deg
        ha_text = 'left'
        text_x = x_start
    else: 
        x_end = lon1 - (lon1 - lon0) * 0.05
        x_start = x_end - bar_width_deg
        ha_text = 'right'
        text_x = x_end

    y_pos = lat0 + (lat1 - lat0) * location_y
    
    ax.plot([x_start, x_end], [y_pos, y_pos], transform=ccrs.PlateCarree(),
            color='black', linewidth=3, zorder=10, clip_on=False)
    
    ax.text(text_x, y_pos - (lat1 - lat0) * 0.05, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha=ha_text, va='top', 
            fontsize=200, clip_on=False)


def plot_yield_and_runoff(file_name, data_dir, fig_path, basin):
    outputs = GetSimSummary(file_name)
    mask = outputs[0]
    data_list = outputs[1:] 
    
    # 1. Load Low Runoff Mask
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        lr_area = ds_lr["Low_Runoff"].where(mask.notnull())

    range_keys = [
        "Yield_kg_ha", "Yield_kg_ha", "Production_ktons",
        "N_Runoff_kg_ha", "N_Runoff_kg_ha", "N_Runoff_kg_ha", "N_Runoff_ktons", "N_Runoff_ktons",
        "P_Runoff_kg_ha", "P_Runoff_kg_ha", "P_Runoff_kg_ha", "P_Runoff_ktons", "P_Runoff_ktons"
    ]
    
    cmaps = ['YlGn']*3 + ['RdPu']*5 + [custom_phosphorus_cmap]*5
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    gdf_boundary = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    
    # Calculate aspect ratio of the basin to force all subplots to be identical
    basin_aspect = (lat_max - lat_min) / (lon_max - lon_min)

    fig, axes = plt.subplots(nrows=13, ncols=1, figsize=(25, 110), 
                             subplot_kw={'projection': ccrs.PlateCarree()})
    
    ims = []
    for i, ax in enumerate(axes):
        vmin, vmax = VAR_RANGES[range_keys[i]]
        
        # FORCE ASPECT: This ensures all maps stay in a perfect vertical line
        ax.set_box_aspect(basin_aspect)
        
        im = data_list[i].plot(ax=ax, transform=ccrs.PlateCarree(), 
                               cmap=cmaps[i], vmin=vmin, vmax=vmax, add_colorbar=False)
        ims.append(im)
        
        # Low Runoff Mask
        lr_area.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=grey_cmap,  add_colorbar=False, zorder=2)
        
        gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=2.0, zorder=3)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.axis('off')

        if i == 12:
            bar_len = 500 if (lon_max - lon_min) > 10 else 100
            add_scale_bar(ax, bar_len, basin)

    # --- Grouped Colorbar Logic with Strict Alignment ---
    cbar_configs = [
        ([0, 1], "Yield [kg/ha]"),
        ([2], "Total Production [ktons]"),
        ([3, 4, 5], "N Runoff [kg/ha]"),
        ([6, 7], "N Runoff [ktons]"),
        ([8, 9, 10], "P Runoff [kg/ha]"),
        ([11, 12], "P Runoff [ktons]")
    ]

    # We determine a fixed 'X' coordinate for all colorbars based on the first map
    # This ensures they are all perfectly vertically aligned
    sample_pos = axes[0].get_position()
    cbar_x_start = sample_pos.x1 + 0.05  # Constant gap from the right edge of the maps

    for indices, label in cbar_configs:
        # Determine top and bottom positions of the group
        pos_top = axes[indices[0]].get_position()
        pos_bot = axes[indices[-1]].get_position()
        
        # Create a new axes for the colorbar at a FIXED X-position
        # [left, bottom, width, height]
        cax = fig.add_axes([cbar_x_start, pos_bot.y0, 0.015, pos_top.y1 - pos_bot.y0])
        
        cbar = fig.colorbar(ims[indices[0]], cax=cax, orientation='vertical', extend='both')
        cbar.set_label(label, fontsize=80, labelpad=80) # Increased labelpad
        cbar.ax.tick_params(labelsize=80, length=20, width=5) # Made ticks bigger to match
        cbar.outline.set_visible(False)

    plt.subplots_adjust(hspace=0.2) 
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

# --- Execution ---
for basin in Studyareas:
    for crop in InputCrops:
        file_name = os.path.join(input_dir, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(file_name):
            print(f"Skipping: {file_name} (Not found)")
            continue
            
        print(f"Processing: {crop} in {basin}...")
        output_file = os.path.join(fig_base_dir, f"{basin}_{crop}_summary.png")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        plot_yield_and_runoff(file_name, data_dir, output_file, basin)

print("Done!")