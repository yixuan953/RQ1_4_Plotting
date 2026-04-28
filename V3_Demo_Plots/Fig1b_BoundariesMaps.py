import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
from matplotlib.colors import ListedColormap, BoundaryNorm

# --- Configuration & Nature Style ---
# Using Liberation Sans as the reliable Linux cluster alternative to Arial
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42  # For editable vector text

Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
InputCrops = ["winterwheat", "mainrice", "secondrice", "soybean", "maize"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
fig_base_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig1_Boundary/1b_Maps"

# --- Define Custom Discrete Colormaps & Norms ---

# Grey Colormap for areas with zero/low runoff
grey_cmap = ListedColormap(["#EBE9E9"])

# 1. Custom Nitrogen (N) Colormap Configuration

n_colors_hex_final = ["#CD2C58", "#E06B80", "#F1998F", "#FFC69D", "#FFE6D4"]

n_boundaries = [0, 5, 15, 20, 25, 50] # 50 is an arbitrary large max value
n_cmap = ListedColormap(n_colors_hex_final[::-1]) # [::-1] Reverse list to go low->high
n_norm = BoundaryNorm(n_boundaries, n_cmap.N)

# 2. Custom Phosphorus (P) Colormap Configuration
p_colors_hex_final = ["#FA812F", "#FFB22C", "#F3C623", "#F9DC85", "#FEF3E2"]

p_boundaries = [0, 0.25, 0.5, 0.75, 1.0, 5.0] # 5.0 is an arbitrary large max value
p_cmap = ListedColormap(p_colors_hex_final[::-1]) # [::-1] Reverse list to go low->high
p_norm = BoundaryNorm(p_boundaries, p_cmap.N)

# Combine configurations (We still create individual plots per unit, but share cmap/norm)
VAR_CONFIG = {
    "N_Runoff_kg_ha": {"cmap": n_cmap, "norm": n_norm, "label": "N Boundary ($kg\ N\ ha^{-1}$)", "ticks": [0, 5, 15, 20, 25, 50], "tick_labels": ['0', '5', '15', '20', '25', '']},
    "N_Runoff_ktons": {"cmap": n_cmap, "norm": n_norm, "label": "N Boundary (ktons)", "ticks": [0, 0.5, 1.0, 1.5, 2.0, 5.0], "tick_labels": ['0', '0.5', '1.0', '1.5', '2.0', '']},
    "P_Runoff_kg_ha": {"cmap": p_cmap, "norm": p_norm, "label": "P Boundary ($kg\ P\ ha^{-1}$)", "ticks": [0, 0.25, 0.5, 0.75, 1.0, 5.0], "tick_labels": ['0', '0.25', '0.5', '0.75', '1.0', '']},
    "P_Runoff_ktons": {"cmap": p_cmap, "norm": p_norm, "label": "P Boundary (ktons)", "ticks": [0, 0.005, 0.01, 0.015, 0.02, 0.50], "tick_labels": ['0', '0.005', '0.01', '0.015', '0.02', '']}
}

# --- Standard Utility Functions (Unchanged from previous update) ---
def GetSimSummary(file_name):
    ds = xr.open_dataset(file_name)
    mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
    
    data_dict = {
        "N_Runoff_kg_ha": (ds["Crit_N_Runoff"] / ds["Total_HA"]) * mask,
        "N_Runoff_ktons": 1e-6 * ds["Crit_N_Runoff"] * mask,
        "P_Runoff_kg_ha": (ds["Crit_P_Runoff"] / ds["Total_HA"]) * mask,
        "P_Runoff_ktons": 1e-6 * ds["Crit_P_Runoff"] * mask
    }
    return mask, data_dict

def add_scale_bar(ax, length_km, basin):
    """Adds a scale bar that scales correctly with dynamic widths."""
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    
    # Position logic
    x_start = lon0 + (lon1 - lon0) * 0.05 if basin == "Yangtze" else lon1 - (lon1 - lon0) * 0.05 - bar_width_deg
    y_pos = lat0 + (lat1 - lat0) * 0.05
    
    ax.plot([x_start, x_start + bar_width_deg], [y_pos, y_pos], 
            transform=ccrs.PlateCarree(), color='black', linewidth=1.5, zorder=10)
    ax.text(x_start + bar_width_deg/2, y_pos + (lat1-lat0)*0.015, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha='center', va='bottom', fontsize=20)

# --- Updated Plotting Function ---
def plot_single_map(data_array, lr_mask, gdf_boundary, var_key, basin, crop, save_path):
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    aspect = (lon_max - lon_min) / (lat_max - lat_min)
    
    # Adjust height to give room for the bottom colorbar
    fixed_height = 4.0 
    dynamic_width = fixed_height * aspect
    
    fig = plt.figure(figsize=(dynamic_width, fixed_height))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    cfg = VAR_CONFIG[var_key]
    
    # 1. Main Data Plot
    im = data_array.plot(ax=ax, cmap=cfg["cmap"], norm=cfg["norm"], add_colorbar=False, zorder=1)
    
    # 2. Low Runoff Mask
    lr_mask.plot(ax=ax, cmap=grey_cmap, add_colorbar=False, zorder=2)
    
    # 3. Basin Boundary
    gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=0.8, zorder=3)
    
    # 4. Scale Bar
    bar_len = 500 if (lon_max - lon_min) > 10 else 100
    add_scale_bar(ax, bar_len, basin)
    
    # 5. Colorbar at the Bottom
    # 'pad' controls distance from map; 'fraction' and 'aspect' control bar size
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.08, fraction=0.046, aspect=30)
    
    # Set custom ticks and labels
    cbar.set_ticks(cfg["ticks"])
    cbar.ax.set_xticklabels(cfg["tick_labels"], fontsize=16)
    cbar.set_label(cfg["label"], fontsize=18, labelpad=10)
    
    # Nature style: remove outer box
    cbar.outline.set_visible(False) 

    ax.set_extent([lon_min, lon_max, lat_min, lat_max])
    ax.axis('off')
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

# --- Execution ---
# Assumes input_dir, data_dir, fig_base_dir defined above are valid
stats_results = []

for basin in Studyareas:
    shp_path = os.path.join(data_dir, "2_shp_StudyArea", basin, f"{basin}.shp")
    if not os.path.exists(shp_path): continue
    gdf_boundary = gpd.read_file(shp_path)
    
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    if not os.path.exists(low_runoff_path): continue
    ds_lr = xr.open_dataset(low_runoff_path)

    for crop in InputCrops:
        file_path = os.path.join(input_dir, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(file_path): continue
        
        print(f"Processing {basin} - {crop}...")
        mask, data_vars = GetSimSummary(file_path)
        lr_mask = ds_lr["Low_Runoff"].where(mask.notnull())

        for var_name, data_array in data_vars.items():
            # Generate Individual Plot (Organized by Basin/Crop)
            out_name = f"{basin}_{crop}_{var_name}.png"
            out_path = os.path.join(fig_base_dir, basin, crop, out_name)
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            
            plot_single_map(data_array, lr_mask, gdf_boundary, var_name, basin, crop, out_path)

            # Statistics (Ignoring 0 and NaN)
            vals = data_array.values.flatten()
            vals = vals[~np.isnan(vals) & (vals > 0)]
            if len(vals) > 0:
                stats_results.append({
                    "Basin": basin, "Crop": crop, "Variable": var_name,
                    "25th": np.percentile(vals, 25), "Median": np.median(vals), "75th": np.percentile(vals, 75)
                })

# Save stats
stats_csv_path = os.path.join(fig_base_dir, "Summary_Statistics.csv")
pd.DataFrame(stats_results).to_csv(stats_csv_path, index=False)
print(f"Done! Maps and statistics have been generated. Stats saved to: {stats_csv_path}")