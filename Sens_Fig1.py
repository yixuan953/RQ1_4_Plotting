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
# Include both mainrice and secondrice explicitly to catch all spatial files
AllCrops = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Sensitivity/N25P10"
os.makedirs(fig_base_dir, exist_ok=True)

N_crit_conc_org = 2.5 
P_crit_conc_org = 0.1

N_crit_conc_sens  = 2.5 # [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
P_crit_conc_sens  = 0.10 # [0.025, 0.05, 0.075, 0.1, 0.125, 0.16]

N_adj_factor = N_crit_conc_sens / N_crit_conc_org
P_adj_factor = P_crit_conc_sens / P_crit_conc_org

# --- 2. Color Schemes & Config ---
grey_cmap = ListedColormap(["#EBE9E9"])

n_colors = ["#D7EBF3", "#FDF2EA", "#FAD2B5", "#F7B296", "#F1998F", "#EB8287", "#E06B80", "#D94C6D", "#CD2C58", "#981E42"]
n_boundaries = [-5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
n_cmap = ListedColormap(n_colors)
n_norm = BoundaryNorm(n_boundaries, n_cmap.N)
n_major_ticks = [0, 10, 20, 30, 40] 

p_colors = ["#D7EBF3", "#FEF3E2", "#F9DC85", "#F6D155", "#F3C623", "#F9C026", "#FFB22C", "#FF992C", "#F17727", "#DF6418"]
p_boundaries = [-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
p_cmap = ListedColormap(p_colors)
p_norm = BoundaryNorm(p_boundaries, p_cmap.N)
p_major_ticks = [0, 0.5, 1.0, 1.5, 2.0]

VAR_CONFIG = {
    "N_Runoff_kg_ha": {
        "cmap": n_cmap, "norm": n_norm, "ticks": n_major_ticks,
        "tick_labels": [str(t) for t in n_major_ticks]
    },
    "P_Runoff_kg_ha": {
        "cmap": p_cmap, "norm": p_norm, "ticks": p_major_ticks,
        "tick_labels": [f"{t:.1f}" if t % 1 != 0 else str(int(t)) for t in p_major_ticks]
    }
}

# --- 3. Helper Functions ---
def add_scale_bar(ax, length_km, basin):
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    
    if basin == "Yangtze":
        x_start_frac = 0.02 
    elif basin == "Rhine":
        bar_frac = bar_width_deg / (lon1 - lon0)
        x_start_frac = 1.25 - bar_frac
    else:
        bar_frac = bar_width_deg / (lon1 - lon0)
        x_start_frac = 0.95 - bar_frac

    x_start_lon = lon0 + (lon1 - lon0) * x_start_frac
    y_pos_lat = lat0 + (lat1 - lat0) * 0.05
    
    ax.plot([x_start_lon, x_start_lon + bar_width_deg], [y_pos_lat, y_pos_lat], 
            transform=ccrs.PlateCarree(), color='black', linewidth=2.0, zorder=10, clip_on=False)
    text_x = x_start_lon + (bar_width_deg / 2)
    ax.text(text_x, y_pos_lat + (lat1-lat0)*0.02, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha='center', va='bottom', fontsize=30)

def plot_total_basin(basin, var_type):
    """Generates a single comprehensive map aggregating all crops for a basin."""
    v_suffix = "N_Runoff_kg_ha" if var_type == 'N' else "P_Runoff_kg_ha"
    sens_adj_factor = N_adj_factor if var_type == 'N' else P_adj_factor
    color_cfg = VAR_CONFIG[v_suffix]
    
    # Setup spatial boundary
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    gdf_boundary = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    aspect = (lon_max - lon_min) / (lat_max - lat_min)

    # Single plot setup
    fig_height = 6.0
    fig, ax = plt.subplots(figsize=(fig_height * aspect, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})

    # Load low runoff mask
    lr_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    ds_lr_all = xr.open_dataset(lr_path)

    # --- AGGREGATION LOGIC ACROSS ALL CROPS ---
    running_total_ha = None
    running_actual_runoff = None
    running_crit_runoff = None
    base_mask = None

    for crop in AllCrops:
        p = os.path.join(input_dir, f"{basin}_{crop}_summary.nc")
        if os.path.exists(p):
            with xr.open_dataset(p) as ds:
                # Save the underlying basin grid layout from the first available file
                if base_mask is None:
                    base_mask = ds["Basin_mask"]

                # Extract absolute values (or calculate absolute loads if variables represent rates)
                # If your files store absolute loads, use them directly:
                actual_crop = ds[f"{var_type}_Runoff"].fillna(0)
                crit_crop = sens_adj_factor * ds[f"Crit_{var_type}_Runoff"].fillna(0)
                ha_crop = ds["Total_HA"].fillna(0)

                if running_total_ha is None:
                    running_total_ha = ha_crop
                    running_actual_runoff = actual_crop
                    running_crit_runoff = crit_crop
                else:
                    running_total_ha += ha_crop
                    running_actual_runoff += actual_crop
                    running_crit_runoff += crit_crop

    if running_total_ha is not None and running_total_ha.max() > 0:
        # Combined rates: Total Load / Total HA
        actual_total = running_actual_runoff / running_total_ha
        crit_total = running_crit_runoff / running_total_ha
        
        # Apply HA threshold filter (> 2500 ha combined across all crops)
        final_mask = base_mask.where(running_total_ha > 2500)
        data_to_plot = (actual_total - crit_total) * final_mask

        # --- PLOTTING ---
        gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=1.0, zorder=3)
        
        im = data_to_plot.plot(ax=ax, cmap=color_cfg["cmap"], norm=color_cfg["norm"], 
                               add_colorbar=False, zorder=1)
        
        # Mask low runoff areas
        current_lr_mask = ds_lr_all["Low_Runoff"].reindex_like(data_to_plot, method='nearest').where(data_to_plot.notnull())
        current_lr_mask.plot(ax=ax, cmap=grey_cmap, add_colorbar=False, zorder=2)

        # Scale Bar
        bar_len = 500 if (lon_max - lon_min) > 10 else 200
        add_scale_bar(ax, bar_len, basin)
        
        # Colorbar
        cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.08, shrink=0.7)
        cbar.set_ticks(color_cfg["ticks"])
        cbar.ax.set_xticklabels(color_cfg["tick_labels"], fontsize=20)
        cbar.outline.set_visible(False)
        
    else:
        ax.text(0.5, 0.5, "No Data Found", transform=ax.transAxes, ha='center', va='center', fontsize=20)

    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.axis('off')

    out_path = os.path.join(fig_base_dir, f"Total_all_crops_exceedance_{basin}_{var_type}.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

# --- Run ---
for basin in Studyareas:
    print(f"Processing Total Basin Maps for: {basin}")
    plot_total_basin(basin, 'N')
    plot_total_basin(basin, 'P')

print("Plotting complete.")