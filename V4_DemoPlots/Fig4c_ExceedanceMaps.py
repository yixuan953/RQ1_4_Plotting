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
import cartopy.geodesic as cgeo
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.gridspec as gridspec

# --- 1. Configuration ---
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
CropGroups = ["winterwheat", "maize", "Rice", "soybean"] 
TargetCrops = ["winterwheat", "maize", "Rice", "soybean"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig4"
os.makedirs(fig_base_dir, exist_ok=True)

# --- 2. Color Schemes & Config ---
# Grey Colormap for areas with zero/low runoff
grey_cmap = ListedColormap(["#EBE9E9"])

# List starts with Blue (Safe), ends with Red (Danger)
n_colors = [
    "#D7EBF3", # Safe/Low
    "#FDF2EA", "#FAD2B5", "#F7B296", "#F1998F", 
    "#EB8287", "#E06B80", "#D94C6D", "#CD2C58", "#981E42" # High
]

# 10 colors, 11 boundaries. 
# We start at -5 so 0 is clearly in the Blue bin.
n_boundaries = [-5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
n_cmap = ListedColormap(n_colors) # REMOVED [::-1]
n_norm = BoundaryNorm(n_boundaries, n_cmap.N)
n_major_ticks = [0, 10, 20, 30, 40] 

# --- 2. Phosphorus (P) ---
p_colors = [
    "#D7EBF3", # Safe/Low
    "#FEF3E2", "#F9DC85", "#F6D155", "#F3C623", 
    "#F9C026", "#FFB22C", "#FF992C", "#F17727", "#DF6418" # High
]

# 10 colors, 11 boundaries
p_boundaries = [-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
p_cmap = ListedColormap(p_colors) # REMOVED [::-1]
p_norm = BoundaryNorm(p_boundaries, p_cmap.N)
p_major_ticks = [0, 0.5, 1.0, 1.5, 2.0]

VAR_CONFIG = {
    "N_Runoff_kg_ha": {
        "cmap": n_cmap, 
        "norm": n_norm, 
        "label": "N Boundary ($kg\ N\ ha^{-1}$)", 
        "ticks": n_major_ticks,
        "tick_labels": [str(t) for t in n_major_ticks]
    },
    "P_Runoff_kg_ha": {
        "cmap": p_cmap, 
        "norm": p_norm, 
        "label": "P Boundary ($kg\ P\ ha^{-1}$)", 
        "ticks": p_major_ticks,
        "tick_labels": [f"{t:.1f}" if t % 1 != 0 else str(int(t)) for t in p_major_ticks]
    }
}


# --- 3. Helper Functions ---

def add_scale_bar(ax, length_km, basin):
    """Adds a scale bar with basin-specific offsets to avoid overlaps."""
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    
    # Calculate how many degrees of longitude equal the requested km at this latitude
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    
    # Position Logic using Axes Fractions (0 to 1)
    # y_pos is 5% from the bottom
    # x_pos varies by basin
    if basin == "Yangtze":
        # Move Yangtze further LEFT (e.g., 2% from the left edge)
        x_start_frac = 0.02 
        ha_align = 'left'
    elif basin == "Rhine":
        # Move Rhine to the bottom LEFT as requested
        bar_frac = bar_width_deg / (lon1 - lon0)
        x_start_frac = 1.25 - bar_frac
        ha_align = 'right'
    else:
        bar_frac = bar_width_deg / (lon1 - lon0)
        x_start_frac = 0.95 - bar_frac
        ha_align = 'center'

    # Convert fraction back to longitude for plotting
    x_start_lon = lon0 + (lon1 - lon0) * x_start_frac
    y_pos_lat = lat0 + (lat1 - lat0) * 0.05
    
    # Plot the line
    ax.plot([x_start_lon, x_start_lon + bar_width_deg], [y_pos_lat, y_pos_lat], 
            transform=ccrs.PlateCarree(), color='black', linewidth=2.0, zorder=10, clip_on=False)
    
    # Plot the text
    text_x = x_start_lon + (bar_width_deg / 2)
    ax.text(text_x, y_pos_lat + (lat1-lat0)*0.02, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha='center', va='bottom', fontsize=30)

def plot_composite_basin(basin, var_type):
    """
    Generates a 1 column x 4 row figure for a specific basin.
    var_type: 'N' or 'P'
    """
    # 1. Setup Color Config
    v_suffix = "N_Runoff_kg_ha" if var_type == 'N' else "P_Runoff_kg_ha"
    color_cfg = VAR_CONFIG[v_suffix]
    
    # 2. Setup Spatial Boundary
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    gdf_boundary = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    aspect = (lon_max - lon_min) / (lat_max - lat_min)

    # 3. Create Figure
    row_height = 5.0
    fig_width = row_height * aspect
    fig = plt.figure(figsize=(fig_width, row_height * len(TargetCrops)))
    gs = gridspec.GridSpec(len(TargetCrops), 1, figure=fig, hspace=0.05)

    # Load shared mask
    lr_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    ds_lr_all = xr.open_dataset(lr_path)

    last_im = None

    for i, group in enumerate(TargetCrops):
        ax = fig.add_subplot(gs[i], projection=ccrs.PlateCarree())
        
        # --- DATA PROCESSING LOGIC ---
        data_to_plot = None
        current_lr_mask = None

        if group == "Rice":
            sc_list = ["mainrice", "secondrice"] if basin == "Yangtze" else ["mainrice"]
            datasets = []
            for sc in sc_list:
                p = os.path.join(input_dir, f"{basin}_{sc}_summary.nc")
                if os.path.exists(p): datasets.append(xr.open_dataset(p))
            
            if datasets:
                master_ds = datasets[0]
                ds_combined = master_ds.fillna(0).copy()
                for v in ["Total_HA", "Crit_N_Runoff", "Crit_P_Runoff"]:
                    ds_combined[v] = sum(d[v].fillna(0) for d in datasets)
                
                # Calculate Exceedance
                crit = (ds_combined[f"Crit_{var_type}_Runoff"]/ds_combined["Total_HA"]).where(ds_combined["Total_HA"] > 0)
                actual = (sum(d[f"{var_type}_Runoff"].fillna(0) for d in datasets) / 
                          ds_combined["Total_HA"]).where(ds_combined["Total_HA"] > 0)
                
                final_mask = master_ds["Basin_mask"].where(ds_combined["Total_HA"] > 2500)
                data_to_plot = (actual - crit) * final_mask
        else:
            p = os.path.join(input_dir, f"{basin}_{group}_summary.nc")
            if os.path.exists(p):
                ds = xr.open_dataset(p)
                crit = (ds[f"Crit_{var_type}_Runoff"]/ds["Total_HA"]).where(ds["Total_HA"] > 0)
                actual = (ds[f"{var_type}_Runoff"]/ds["Total_HA"]).where(ds["Total_HA"] > 0) # This is kg/ha in your summary files
                final_mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500)
                data_to_plot = (actual - crit) * final_mask

        # --- PLOTTING ---
        gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=0.8, zorder=3)
        
        if data_to_plot is not None:
            im = data_to_plot.plot(ax=ax, cmap=color_cfg["cmap"], norm=color_cfg["norm"], 
                                   add_colorbar=False, zorder=1)
            
            # Mask low runoff
            current_lr_mask = ds_lr_all["Low_Runoff"].reindex_like(data_to_plot, method='nearest').where(data_to_plot.notnull())
            current_lr_mask.plot(ax=ax, cmap=grey_cmap, add_colorbar=False, zorder=2)
            last_im = im


        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.axis('off')

        # Scale bar on bottom plot (i=3)
        if i == len(TargetCrops) - 4:
            bar_len = 500 if (lon_max - lon_min) > 10 else 200
            add_scale_bar(ax, bar_len, basin)

    # 4. Shared Colorbar
    if last_im:
        cbar_ax = fig.add_axes([0.2, 0.07, 0.6, 0.015])
        cbar = fig.colorbar(last_im, cax=cbar_ax, orientation='horizontal')
        cbar.set_ticks(color_cfg["ticks"])
        cbar.ax.set_xticklabels(color_cfg["tick_labels"], fontsize=20)
        cbar.set_label(f"{var_type} Exceedance ($kg\ {var_type}\ ha^{{-1}}$)", fontsize=22, labelpad=10)
        cbar.outline.set_visible(False)

    # Save
    out_path = os.path.join(fig_base_dir, f"Boundary_exceedance_maps_{basin}_{var_type}.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

# --- Run ---
for basin in Studyareas:
    print(f"Processing Basin: {basin}")
    plot_composite_basin(basin, 'N')
    plot_composite_basin(basin, 'P')

print("Plotting complete.")