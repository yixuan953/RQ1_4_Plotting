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
import cartopy.geodesic as cgeo
from matplotlib.colors import ListedColormap, BoundaryNorm

# --- Configuration ---
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']

Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
DisplayCrops = ["Wheat", "Maize", "Rice", "Soybean"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig2_Sim_Yield_Loss/2b_Maps"

# --- Colormaps & Norms ---
grey_cmap = ListedColormap(["#EBE9E9"])

# 1. Production (Green)
prod_colors_hex = [
    "#FDFBD3", # 0-100 (Cream)
    "#F8F2B2", # 100-200
    "#F0E491", # 200-300 (Original Light Yellow)
    "#D6D67A", # 300-400
    "#BBC863", # 400-500 (Original Lime/Olive)
    "#8FAB5D", # 500-600
    "#658C58", # 600-700 (Original Mid Green)
    "#4B7A53", # 700-800
    "#31694E"  # 800-1000 (Original Deep Green)
]
prod_cmap = ListedColormap(prod_colors_hex)
prod_bounds = [0, 100, 200, 300, 400, 500, 600, 700, 800, 1000] 
prod_norm = BoundaryNorm(prod_bounds, prod_cmap.N)

# 2. Nitrogen (Pink/Red)
n_colors = ["#FFE6D4", "#FFC69D", "#F1998F", "#E06B80", "#CD2C58"]
n_bounds = [0, 10, 20, 30, 40, 100]
n_cmap = ListedColormap(n_colors)
n_norm = BoundaryNorm(n_bounds, n_cmap.N)

# 3. Phosphorus (Orange/Yellow)
p_colors = ["#FEF3E2", "#F9DC85", "#F3C623", "#FFB22C", "#FA812F"]
p_bounds = [0, 0.5, 1.0, 1.5, 2.0, 5.0]
p_cmap = ListedColormap(p_colors)
p_norm = BoundaryNorm(p_bounds, p_cmap.N)

VAR_CONFIG = {
    "Production_ktons": {
        "cmap": prod_cmap, 
        "norm": prod_norm, 
        "label": "Production (ktons)", 
        "ticks": prod_bounds, 
        "tick_labels": ['0', '', '0.2', '', '0.4', '', '0.6', '', '0.8',''] 
    },

    "N_Runoff_kg_ha": {
        "cmap": n_cmap, "norm": n_norm, "label": "N Runoff ($kg\ N\ ha^{-1}$)", 
        "ticks": n_bounds, "tick_labels": ['0', '10', '20', '30', '40', '']
    },
    "P_Runoff_kg_ha": {
        "cmap": p_cmap, "norm": p_norm, "label": "P Runoff ($kg\ P\ ha^{-1}$)", 
        "ticks": p_bounds, "tick_labels": ['0', '0.5', '1.0', '1.5', '2.0', '']
    }
}

# --- Helper Functions ---

def get_crop_files(basin, display_crop):
    if display_crop == "Rice":
        return ["mainrice", "secondrice"] if basin == "Yangtze" else ["mainrice"]
    if display_crop == "Wheat": return ["winterwheat"]
    return [display_crop.lower()]

def aggregate_data(basin, suffixes):
    t_prod, t_n_mass, t_p_mass, t_ha, master_mask = None, None, None, None, None
    
    for s in suffixes:
        path = os.path.join(input_dir, f"{basin}_{s}_summary.nc")
        if not os.path.exists(path): 
            continue
            
        with xr.open_dataset(path) as ds:
            # 1. Create the specific mask for this file
            mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
            
            # 2. Update master_mask to be the union of all crop areas
            if master_mask is None: 
                master_mask = mask
            else:
                # This ensures master_mask covers everywhere ANY crop is grown
                master_mask = master_mask.fillna(mask)
            
            # 3. Calculate variables for THIS crop
            prod = 1e-6 * (ds["Avg_Yield_Irrigated"].fillna(0) * ds["Irrigated_HA"].fillna(0) + 
                           ds["Avg_Yield_Rainfed"].fillna(0) * ds["Rainfed_HA"].fillna(0))
            n_mass = ds["N_Runoff"].fillna(0)
            p_mass = ds["P_Runoff"].fillna(0)
            ha = ds["Total_HA"].fillna(0)

            # 4. Aggregate using fillna(0) to prevent NaN "poisoning"
            if t_prod is None:
                t_prod, t_n_mass, t_p_mass, t_ha = prod, n_mass, p_mass, ha
            else:
                t_prod = t_prod.fillna(0) + prod.fillna(0)
                t_n_mass = t_n_mass.fillna(0) + n_mass.fillna(0)
                t_p_mass = t_p_mass.fillna(0) + p_mass.fillna(0)
                t_ha = t_ha.fillna(0) + ha.fillna(0)

    if t_ha is None:
        return None, None

    # 5. Final Step: Apply the master_mask so we only show the basin area
    # and ensure we don't divide by zero
    results = {
        "Production_ktons": t_prod.where(master_mask.notnull()),
        "N_Runoff_kg_ha": (t_n_mass / t_ha).where(t_ha > 0).where(master_mask.notnull()),
        "P_Runoff_kg_ha": (t_p_mass / t_ha).where(t_ha > 0).where(master_mask.notnull())
    }

    return results, master_mask

def add_scale_bar(ax, length_km, basin):
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    x_start = lon0 + (lon1 - lon0) * 0.05 if basin == "Yangtze" else lon1 - (lon1 - lon0) * 0.05 - bar_width_deg
    y_pos = lat0 + (lat1 - lat0) * 0.05
    ax.plot([x_start, x_start + bar_width_deg], [y_pos, y_pos], transform=ccrs.PlateCarree(), color='black', linewidth=1.5, zorder=10)
    ax.text(x_start + bar_width_deg/2, y_pos + (lat1-lat0)*0.015, f'{length_km} km', transform=ccrs.PlateCarree(), ha='center', va='bottom', fontsize=18)

def plot_single_map(data_array, lr_mask, gdf_boundary, var_key, basin, crop, save_path):
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    aspect = (lon_max - lon_min) / (lat_max - lat_min)
    fixed_height = 4.0 
    fig = plt.figure(figsize=(fixed_height * aspect, fixed_height))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    cfg = VAR_CONFIG[var_key]
    im = data_array.plot(ax=ax, cmap=cfg["cmap"], norm=cfg["norm"], add_colorbar=False, zorder=1)
    lr_mask.plot(ax=ax, cmap=grey_cmap, add_colorbar=False, zorder=2)
    gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=0.8, zorder=3)
    
    add_scale_bar(ax, (500 if (lon_max - lon_min) > 10 else 100), basin)
    
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.08, fraction=0.046, aspect=30)
    cbar.set_ticks(cfg["ticks"])
    cbar.ax.set_xticklabels(cfg["tick_labels"], fontsize=16)
    cbar.set_label(cfg["label"], fontsize=18)
    cbar.outline.set_visible(False) 

    ax.set_extent([lon_min, lon_max, lat_min, lat_max])
    ax.axis('off')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

# --- Execution ---
for basin in Studyareas:
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    if not os.path.exists(shp_path): continue
    gdf_boundary = gpd.read_file(shp_path)
    
    lr_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    ds_lr = xr.open_dataset(lr_path)

    for crop in DisplayCrops:
        suffixes = get_crop_files(basin, crop)
        data_dict, mask = aggregate_data(basin, suffixes)
        if mask is None: continue
        
        lr_mask = ds_lr["Low_Runoff"].where(mask.notnull())

        for var_name, data_array in data_dict.items():
            out_dir = os.path.join(fig_base_dir, basin, crop)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"{basin}_{crop}_{var_name}.png")
            
            plot_single_map(data_array, lr_mask, gdf_boundary, var_name, basin, crop, out_path)
            print(f"Saved: {basin} {crop} {var_name}")