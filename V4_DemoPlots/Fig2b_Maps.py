import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm

# --- 1. Configuration ---
STUDY_AREAS = ["LaPlata", "Indus", "Yangtze", "Rhine"]
TARGET_CROPS = ["Wheat", "Maize", "Rice", "Soybean"]

INPUT_DIR = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
DATA_DIR = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
FIG_OUT = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig2"
os.makedirs(FIG_OUT, exist_ok=True)

# --- 2. Custom Color Setup ---
grey_cmap = ListedColormap(["#ece9e9"])

# Nitrogen (Blue to Red) - 10 colors
n_colors = [ "#FDF2EA", "#FAD2B5", "#F7B296", "#F1998F", "#EB8287", "#E06B80", "#D94C6D", "#CD2C58", "#981E42"]
n_bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
n_cmap = ListedColormap(n_colors)
n_norm = BoundaryNorm(n_bounds, n_cmap.N)

# Phosphorus (Blue to Orange) - 10 colors
p_colors = ["#FEF3E2", "#F9DC85", "#F6D155", "#F3C623", "#F9C026", "#FFB22C", "#FF992C", "#F17727", "#DF6418"]
p_bounds = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
p_cmap = ListedColormap(p_colors)
p_norm = BoundaryNorm(p_bounds, p_cmap.N)

VAR_CONFIG = {
    "N_boundary": {
        "cmap": n_cmap, "norm": n_norm, 
        "label": "N Runoff (kg N/ha)", "ticks": [0, 10, 20, 30, 40]
    },
    "P_boundary": {
        "cmap": p_cmap, "norm": p_norm, 
        "label": "P Runoff (kg P/ha)", "ticks": [0, 0.5, 1.0, 1.5, 2.0]
    }
}

def get_crop_files(basin, crop):
    if crop == "Rice":
        return ["mainrice"]
    if crop == "Wheat": return ["winterwheat"]
    return [crop.lower()]

def aggregate_data(basin, crop_list):
    t_crit_n, t_crit_p, N_boundary, P_boundary, t_ha, master_mask = None, None, None, None, None, None
    for c in crop_list:
        path = os.path.join(INPUT_DIR, f"{basin}_{c}_summary.nc")
        if not os.path.exists(path): continue
        with xr.open_dataset(path) as ds:
            mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
            
            t_crit_n = ds["Crit_N_Runoff"]
            t_crit_p = ds["Crit_P_Runoff"]
            t_ha = ds["Total_HA"]

            N_boundary = mask* t_crit_n/t_ha
            P_boundary = mask* t_crit_p/t_ha

            if N_boundary is None:
                return None

            # Return a dictionary instead of a tuple
            return {
                "N_boundary": N_boundary,
                "P_boundary": P_boundary,
                "t_ha": t_ha,
                "mask": mask
            }


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


for basin in STUDY_AREAS:
    print(f"Processing Basin: {basin}")
    
    # 1. Load Basin files
    shp_path = os.path.join(DATA_DIR, "2_shp_StudyArea", basin, f"{basin}.shp")
    if not os.path.exists(shp_path): continue
    gdf = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf.total_bounds
    aspect = (lon_max - lon_min) / (lat_max - lat_min)

    lr_path = os.path.join(DATA_DIR, "2_StudyArea", basin, "low_runoff_mask.nc")
    ds_lr = xr.open_dataset(lr_path)

    # 2. Iterate through each Variable (One 4-row figure per variable)
    for var_key in ["N_boundary", "P_boundary"]:
        cfg = VAR_CONFIG[var_key]
        
        # Setup Figure: 1 column, 4 rows
        fig = plt.figure(figsize=(6 * aspect, 18)) # Adjust size for verticality
        gs = gridspec.GridSpec(4, 1, figure=fig, hspace=0.05)
        
        last_im = None

        for i, crop in enumerate(TARGET_CROPS):
            ax = fig.add_subplot(gs[i], projection=ccrs.PlateCarree())
            data_dict = aggregate_data(basin, get_crop_files(basin, crop))
            
            if data_dict is not None:
                val = data_dict[var_key]
                
                # Plot Data
                im = val.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cfg["cmap"], norm=cfg["norm"], add_colorbar=False, zorder=1)
                
                # Plot Grey Mask (Low Runoff)
                lr_layer = ds_lr["Low_Runoff"].reindex_like(val, method='nearest')
                lr_layer.where((lr_layer == 1) & (val.notnull())).plot(ax=ax, transform=ccrs.PlateCarree(), cmap=grey_cmap, add_colorbar=False, zorder=2)
                
                last_im = im
            
            # Plot Boundary
            gdf.boundary.plot(ax=ax, color='black', linewidth=0.8, zorder=3)
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax.axis('off')

            # Scale bar on bottom-most plot
            if i == 0:
                s_len = 500 if (lon_max-lon_min) > 10 else 200
                add_scale_bar(ax, s_len, basin)

        # 3. Add Shared Colorbar for the Basin-Variable set
        if last_im:
            cbar_ax = fig.add_axes([0.25, 0.08, 0.5, 0.015]) # [left, bottom, width, height]
            cbar = fig.colorbar(last_im, cax=cbar_ax, orientation='horizontal')
            cbar.set_label(cfg["label"], fontsize=30)
            cbar.set_ticks(cfg["ticks"])
            cbar.ax.tick_params(labelsize=25)
            cbar.outline.set_visible(False)

        # 4. Save the Final Stacked Map
        save_path = os.path.join(FIG_OUT, f"{basin}_{var_key}_Maps.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: {save_path}")