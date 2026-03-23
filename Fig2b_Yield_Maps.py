# This code is used to plot: 
# 1 - Rainfed yield [kg/ha]
# 2 - Irrigated yield [kg/ha]
# 3 - Total crop production [kg/ha]

# 4 - Rainfed N runoff [kg/ha]
# 5 - Irrigated N runoff [kg/ha]
# 6 - Total N runoff [ktons]

# 7 - Rainfed P runoff [kg/ha]
# 8 - Irrigated P runoff [kg/ha]
# 9 - Total P runoff [ktons]


import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import cartopy.geodesic as cgeo

# --- Configuration ---
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
InputCrops = ["winterwheat", "mainrice", "secondrice", "soybean", "maize"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig2_Sim_Yield_Loss/2b_Maps" 

DisplayCrops = ["Wheat", "Maize", "Rice", "Soybean"]

# Colormaps
p_colors = ["#fffaf3ff", "#e49503ff"]
custom_phosphorus_cmap = LinearSegmentedColormap.from_list("custom_p", p_colors, N=256)
grey_cmap = ListedColormap(["#EBE9E9"])

# Ranges
VAR_RANGES = {
    "Production_ktons": (0, 500),  
    "N_Runoff_kg_ha":   (0, 40),
    "P_Runoff_kg_ha":   (0, 1.0)
}

def get_crop_files(basin, display_crop):
    if display_crop == "Rice" and basin == "Yangtze":
        return ["mainrice", "secondrice"]
    elif display_crop == "Rice":
        return ["mainrice"]
    elif display_crop == "Wheat":
        return ["winterwheat"]
    else:
        return [display_crop.lower()]

def add_scale_bar(ax, length_km, basin, location_y=-0.08):
    """Adds a scale bar to the map."""
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    
    # Calculate degrees per km at the current latitude
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    
    # Positioning
    if basin == "Yangtze":
        x_start = lon0 + (lon1 - lon0) * 0.05
    else: 
        x_end = lon1 - (lon1 - lon0) * 0.05
        x_start = x_end - bar_width_deg

    y_pos = lat0 + (lat1 - lat0) * location_y
    
    ax.plot([x_start, x_start + bar_width_deg], [y_pos, y_pos], 
            transform=ccrs.PlateCarree(), color='black', linewidth=6, zorder=10, clip_on=False)
    
    ax.text(x_start, y_pos - (lat1 - lat0) * 0.02, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha='left', va='top', 
            fontsize=35, clip_on=False)

def aggregate_crop_data(basin, file_suffixes):
    t_prod, t_n_mass, t_p_mass, t_ha, master_mask = None, None, None, None, None

    for suffix in file_suffixes:
        path = os.path.join(input_dir, f"{basin}_{suffix}_summary.nc")
        if not os.path.exists(path): continue
            
        with xr.open_dataset(path) as ds:
            mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
            if master_mask is None: master_mask = mask
            
            # Production
            prod = 1e-6 * (ds["Avg_Yield_Irrigated"].fillna(0) * ds["Irrigated_HA"].fillna(0) + 
                           ds["Avg_Yield_Rainfed"].fillna(0) * ds["Rainfed_HA"].fillna(0)) * mask
            
            # Runoff logic
            area = ds["Total_HA"] * mask
            n_mass = ds["N_Runoff"] 
            p_mass = ds["P_Runoff"] 

            if t_prod is None:
                t_prod, t_n_mass, t_p_mass, t_ha = prod, n_mass, p_mass, area
            else:
                t_prod = t_prod.fillna(0) + prod.fillna(0)
                t_n_mass = t_n_mass.fillna(0) + n_mass.fillna(0)
                t_p_mass = t_p_mass.fillna(0) + p_mass.fillna(0)
                t_ha = t_ha.fillna(0) + area.fillna(0)

    avg_n = (t_n_mass / t_ha).where(t_ha > 0) * master_mask
    avg_p = (t_p_mass / t_ha).where(t_ha > 0) * master_mask
    t_prod = t_prod * master_mask
    
    return t_prod, avg_n, avg_p, master_mask

def plot_crop_3panel(basin, crop_name, data_package, output_path):
    prod, n_run, p_run, mask = data_package
    
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        lr_area = ds_lr["Low_Runoff"].where(mask.notnull())

    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    gdf_boundary = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    aspect = (lat_max - lat_min) / (lon_max - lon_min)

    fig, axes = plt.subplots(3, 1, figsize=(10, 25), subplot_kw={'projection': ccrs.PlateCarree()})
    
    layers = [prod, n_run, p_run]
    keys = ["Production_ktons", "N_Runoff_kg_ha", "P_Runoff_kg_ha"]
    cmaps = ["YlGn", "RdPu", custom_phosphorus_cmap]

    for i, ax in enumerate(axes):
        vmin, vmax = VAR_RANGES[keys[i]]
        ax.set_box_aspect(aspect)
        
        im = layers[i].plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmaps[i],
                            vmin=vmin, vmax=vmax, add_colorbar=False)
        
        lr_area.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=grey_cmap, add_colorbar=False, zorder=2)
        gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=2.5, zorder=3)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.axis('off')

        # Colorbar
        pos = ax.get_position()
        cax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.03, pos.height])
        cbar = fig.colorbar(im, cax=cax, extend='both')
        cbar.ax.tick_params(labelsize=35)

        # SCALE BAR: Added to the first subplot only
        if i == 0:
            span = lon_max - lon_min
            # Adaptive length based on basin size
            bar_len = 500 if span > 10 else 100 if span > 2 else 50
            add_scale_bar(ax, bar_len, basin)

    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()

# --- Execution ---
for basin in Studyareas:
    for crop in DisplayCrops:
        suffixes = get_crop_files(basin, crop)
        if not os.path.exists(os.path.join(input_dir, f"{basin}_{suffixes[0]}_summary.nc")):
            continue
            
        print(f"Plotting 3-panel maps for {crop} in {basin}...")
        data = aggregate_crop_data(basin, suffixes)
        
        out_name = f"{basin}_{crop}_3Panel_Summary.png"
        plot_crop_3panel(basin, crop, data, os.path.join(fig_base_dir, out_name))