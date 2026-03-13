# This code is used to plot: 
# 1 - Boundary for cropland N runoff to surface water by crop [ktons]
# 2 - Boundary for cropland N runoff to surface water by crop [kg/ha]
# 3 - Boundary for cropland P runoff to surface water by crop [ktons]
# 4 - Boundary for cropland P runoff to surface water by crop [kg/ha]

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
InputCrops =  ["winterwheat", "mainrice", "secondrice", "soybean", "maize"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig1_Boundary/1b_Maps"

# Custom Phosphorus Colormap
p_colors = ["#fffaf3ff", "#ef9b00ff"]
custom_phosphorus_cmap = LinearSegmentedColormap.from_list("custom_p", p_colors, N=256)

# Grey Colormap for areas with low runoff
grey_cmap = ListedColormap(["#EBE9E9"])

VAR_RANGES = {
    "N_Runoff_kg_ha":   (0, 20),
    "N_Runoff_ktons":   (0, 2),
    "P_Runoff_kg_ha":   (0, 0.2),  
    "P_Runoff_ktons":   (0, 0.02)

}


def GetSimSummary(file_name):
    ds = xr.open_dataset(file_name)
    # Threshold masking
    mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
    
    # N Runoff logic (Note: Fix indices/names if they differ in your .nc)
    Crit_N_kgperha     = (ds["Crit_N_Runoff"] / ds["Total_HA"]) * mask
    Crit_N   = 1e-6 * ds["Crit_N_Runoff"] * mask

    # P Runoff logic
    Crit_P_kgperha     = (ds["Crit_P_Runoff"] / ds["Total_HA"]) * mask
    Crit_P   = 1e-6 * ds["Crit_P_Runoff"] * mask
    
    return (mask, Crit_N_kgperha, Crit_N, Crit_P_kgperha, Crit_P)

def calculate_stats(data_array):
    """Calculates median and percentiles for an xarray DataArray, 
    excluding NaNs AND zeros."""
    
    # Flatten the 2D spatial data into a 1D array
    vals = data_array.values.flatten()
    
    # Filter: Keep only values that are NOT NaN and NOT 0
    # Using np.isclose or > 0 is safer for floating point data
    vals = vals[~np.isnan(vals) & (vals > 0)]
    
    # Handle cases where the entire map might be zero or NaN
    if len(vals) == 0:
        return np.nan, np.nan, np.nan
    
    # Calculate statistics on the remaining non-zero data
    q1 = np.percentile(vals, 25)
    median = np.percentile(vals, 50)
    q3 = np.percentile(vals, 75)
    
    return q1, median, q3

def add_scale_bar(ax, length_km, basin, location_y=-0.15):

    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    
    if basin in ["Yangtze"]:
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
            color='black', linewidth=5, zorder=10, clip_on=False)
    
    ax.text(text_x, y_pos - (lat1 - lat0) * 0.05, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha=ha_text, va='top', 
            fontsize=70, clip_on=False)


def plot_boundaries(file_name, data_dir, fig_path, basin):
    outputs = GetSimSummary(file_name)
    mask = outputs[0]
    data_list = outputs[1:] 
    
    # Load Low Runoff Mask
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        lr_area = ds_lr["Low_Runoff"].where(mask.notnull())

    range_keys = ["N_Runoff_kg_ha", "N_Runoff_ktons",  "P_Runoff_kg_ha", "P_Runoff_ktons"]
    
    cmaps = ['RdPu']*2 + [custom_phosphorus_cmap]*2
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    gdf_boundary = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    
    # Calculate aspect ratio of the basin to force all subplots to be identical
    basin_aspect = (lat_max - lat_min) / (lon_max - lon_min)

    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(30, 50),
                             subplot_kw={'projection': ccrs.PlateCarree()})
    
    for i, ax in enumerate(axes):
        vmin, vmax = VAR_RANGES[range_keys[i]]
        ax.set_box_aspect(basin_aspect)
        
        # 1. Plot the main data
        im = data_list[i].plot(ax=ax, transform=ccrs.PlateCarree(), 
                               cmap=cmaps[i], vmin=vmin, vmax=vmax, add_colorbar=False)
        
        # 2. Add masks and boundaries
        lr_area.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=grey_cmap, add_colorbar=False, zorder=2)
        gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=2.0, zorder=3)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.axis('off')

        # Get position of the current axis to place the colorbar exactly next to it
        pos = ax.get_position()
        
        # [left, bottom, width, height]
        # We place it to the right (pos.x1 + offset), same bottom/height as the map
        cbar_width = 0.04  # Significantly wider
        cbar_offset = 0.02 # Gap between map and bar
        cax = fig.add_axes([pos.x1 + cbar_offset, pos.y0, cbar_width, pos.height])
        
        cbar = fig.colorbar(im, cax=cax, orientation='vertical', extend='both')
        
        # Use the specific label for this index from your range_keys or a mapping
        cbar.ax.tick_params(labelsize=85, length=25, width=6)
        cbar.outline.set_visible(True) # Usually looks better for individual bars

        if i == 0:  
            bar_len = 500 if (lon_max - lon_min) > 10 else 100
            add_scale_bar(ax, bar_len, basin)

    # hspace controls the vertical gap between the maps
    plt.subplots_adjust(hspace=0.3) 
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


stats_results = []
# --- Execution ---
for basin in Studyareas:
    for crop in InputCrops:
        file_name = os.path.join(input_dir, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(file_name):
            print(f"Skipping: {file_name} (Not found)")
            continue
            
        print(f"Processing: {crop} in {basin}...")
        output_file = os.path.join(fig_base_dir, f"{basin}_{crop}_boundaries.png")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Plotting function 
        plot_boundaries(file_name, data_dir, output_file, basin)

        # Output summary statistics
        outputs = GetSimSummary(file_name)
        # indices based on your return (mask, N_kg_ha, N_ktons, P_kg_ha, P_ktons)
        data_vars = {
            "N_Runoff_kg_ha": outputs[1],
            "N_Runoff_ktons": outputs[2],
            "P_Runoff_kg_ha": outputs[3],
            "P_Runoff_ktons": outputs[4]
        }
        for var_name, data_array in data_vars.items():
            q1, median, q3 = calculate_stats(data_array)
            
            stats_results.append({
                "Basin": basin,
                "Crop": crop,
                "Variable": var_name,
                "25th_Percentile": q1,
                "Median": median,
                "75th_Percentile": q3
            })


# --- 3. Create DataFrame and Save ---
df_stats = pd.DataFrame(stats_results)

# Display the first few rows
print(df_stats.head())

# Save to CSV for Nature Food Supplementary Table
stats_csv_path = os.path.join(fig_base_dir, "Summary_Statistics.csv")
df_stats.to_csv(stats_csv_path, index=False)
print(f"Statistics saved to {stats_csv_path}")
print("Done!")