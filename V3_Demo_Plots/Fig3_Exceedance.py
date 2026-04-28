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

# --- 1. Configuration ---
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
CropGroups = ["winterwheat", "maize", "Rice", "soybean"] 

# input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig4_Exceedance"
os.makedirs(fig_base_dir, exist_ok=True)

# --- 2. Color Schemes & Config ---
grey_cmap = ListedColormap(["#EBE9E9"])

n_colors = ["#D9F3B9", "#FFE6D4", "#FFC69D", "#F1998F", "#E06B80", "#CD2C58"]
n_bounds = [-100, 0, 10, 20, 30, 40, 100]
n_cmap, n_norm = ListedColormap(n_colors), BoundaryNorm(n_bounds, len(n_colors))

p_colors = ['#D9F3B9', "#FEF3E2", "#F9DC85", "#F3C623", "#FFB22C", "#FA812F"]
p_bounds = [-5.0, 0, 0.5, 1.0, 1.5, 2.0, 5.0]
p_cmap, p_norm = ListedColormap(p_colors), BoundaryNorm(p_bounds, len(p_colors))

VAR_CONFIG = {
    "N_Exceedance": {"cmap": n_cmap, "norm": n_norm, "ticks": n_bounds, "label": "Nitrogen Exceedance [kg N/ha]"},
    "P_Exceedance": {"cmap": p_cmap, "norm": p_norm, "ticks": p_bounds, "label": "Phosphorus Exceedance [kg P/ha]"}
}

# --- 3. Helper Functions ---

def add_scale_bar(ax, length_km, location_y=-0.08):
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    x_start = lon0 + (lon1 - lon0) * 0.05
    y_pos = lat0 + (lat1 - lat0) * location_y
    ax.plot([x_start, x_start + bar_width_deg], [y_pos, y_pos], transform=ccrs.PlateCarree(),
            color='black', linewidth=3, zorder=10, clip_on=False)
    ax.text(x_start, y_pos - (lat1 - lat0) * 0.02, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha='left', va='top', fontsize=12, clip_on=False)

def plot_single_map(data_array, lr_mask, gdf_boundary, var_key, save_path):
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    aspect = (lon_max - lon_min) / (lat_max - lat_min)
    
    fig = plt.figure(figsize=(6.0 * aspect, 6.0))
    ax = plt.axes(projection=ccrs.PlateCarree())
    cfg = VAR_CONFIG[var_key]
    
    im = data_array.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cfg["cmap"], 
                         norm=cfg["norm"], add_colorbar=False, zorder=1)
    
    lr_mask.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=grey_cmap, add_colorbar=False, zorder=2)
    gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=1.2, zorder=3)
    
    # add_scale_bar(ax, (500 if (lon_max - lon_min) > 10 else 100))
    
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.12, fraction=0.05, aspect=30)
    cbar.set_ticks(cfg["ticks"])
    cbar.ax.set_xticklabels([str(b) for b in cfg["ticks"]], fontsize=20)
    cbar.set_label(cfg["label"], fontsize=14, weight='bold')
    cbar.outline.set_visible(False) 

    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.axis('off')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

# --- 4. Main Processing Loop ---

for basin in Studyareas:
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    if not os.path.exists(shp_path): continue
    gdf_boundary = gpd.read_file(shp_path)
    lr_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    
    for group in CropGroups:
        if group == "Rice":
            # Rice logic: Sum both but use mainrice for the mask/grid
            sc_list = ["mainrice", "secondrice"] if basin == "Yangtze" else ["mainrice"]
            datasets = []
            for sc in sc_list:
                p = os.path.join(input_dir, f"{basin}_{sc}_summary.nc")
                if os.path.exists(p): datasets.append(xr.open_dataset(p))
            
            if not datasets: continue
            
            # Use mainrice as the coordinate and mask master
            master_ds = datasets[0] 
            ds_combined = master_ds.fillna(0).copy()
            
            # Combine Mass/Area
            for v in ["Total_HA", "Crit_N_Runoff", "Crit_P_Runoff"]:
                ds_combined[v] = sum(d[v].fillna(0) for d in datasets)
            for v_base in ["N_Runoff", "P_Runoff"]:
                total_mass = sum((d[v_base] * d["Total_HA"]).fillna(0) for d in datasets)
                ds_combined[v_base] = total_mass / ds_combined["Total_HA"].where(ds_combined["Total_HA"] > 0)
            
            # Apply strict spatial filters: Basin + Area > 2500
            # For rice, we use the mainrice mask as requested
            final_mask = master_ds["Basin_mask"].where(ds_combined["Total_HA"] > 2500)
            
        else:
            # Single crop logic
            p = os.path.join(input_dir, f"{basin}_{group}_summary.nc")
            if not os.path.exists(p): continue
            ds_combined = xr.open_dataset(p)
            final_mask = ds_combined["Basin_mask"].where(ds_combined["Total_HA"] > 2500)

        # Calculate Rates and Exceedance with the final mask
        crit_n = (ds_combined["Crit_N_Runoff"] / ds_combined["Total_HA"]).where(ds_combined["Total_HA"] > 0)
        crit_p = (ds_combined["Crit_P_Runoff"] / ds_combined["Total_HA"]).where(ds_combined["Total_HA"] > 0)
        n_avg = (ds_combined["N_Runoff"] / ds_combined["Total_HA"]).where(ds_combined["Total_HA"] > 0)
        p_avg = (ds_combined["P_Runoff"] / ds_combined["Total_HA"]).where(ds_combined["Total_HA"] > 0)

        n_ex = (n_avg - crit_n) * final_mask
        p_ex = (p_avg - crit_p) * final_mask

        # Align and mask Low Runoff
        with xr.open_dataset(lr_path) as ds_lr:
            lr_mask = ds_lr["Low_Runoff"].reindex_like(final_mask, method='nearest').where(final_mask.notnull())

        # Save Maps
        plot_single_map(n_ex, lr_mask, gdf_boundary, "N_Exceedance", os.path.join(fig_base_dir, f"{basin}_{group}_N_Exceedance.png"))
        plot_single_map(p_ex, lr_mask, gdf_boundary, "P_Exceedance", os.path.join(fig_base_dir, f"{basin}_{group}_P_Exceedance.png"))

print("Plotting complete.")