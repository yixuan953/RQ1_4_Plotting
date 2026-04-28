import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
import geopandas as gpd
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm

# --- 1. Configuration ---
STUDY_AREAS = ["LaPlata", "Indus", "Yangtze", "Rhine"]
TARGET_CROPS = ["Wheat", "Maize", "Rice", "Soybean"]
CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}

DATA_DIR = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
BASE_DIR = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
SCEN_DIR = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
FIG_OUT  = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V4_Plots/Fig5"
os.makedirs(FIG_OUT, exist_ok=True)

# --- 2. Custom Color Setup ---
prod_colors = ["#AA2B1D", "#D25D38", "#E89154", "#F7C772", "#E3D160", "#8AA624"]
prod_bounds = [-100, -40, -30, -20, -10, 0, 100]
prod_cmap = mcolors.ListedColormap(prod_colors)
prod_norm = mcolors.BoundaryNorm(prod_bounds, prod_cmap.N)
grey_cmap = mcolors.ListedColormap(["#ece9e9"]) 

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

def get_aggregated_production(folder, basin, crop_list):
    t_prod = None
    for crop in crop_list:
        path = os.path.join(folder, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(path): continue
        with xr.open_dataset(path) as ds:
            mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
            prod = (ds["Avg_Yield_Irrigated"].fillna(0) * ds["Irrigated_HA"].fillna(0) + 
                    ds["Avg_Yield_Rainfed"].fillna(0) * ds["Rainfed_HA"].fillna(0)) * mask
            t_prod = prod if t_prod is None else t_prod.fillna(0) + prod.fillna(0)
    return t_prod

# --- 4. Main Plotting Function ---

def plot_basin_production(basin):
    # 1. Setup Spatial Boundary
    shp_path = os.path.join(DATA_DIR, "2_shp_StudyArea", basin, f"{basin}.shp")
    gdf_boundary = gpd.read_file(shp_path)
    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds
    aspect = (lon_max - lon_min) / (lat_max - lat_min)

    # 2. Setup GridSpec Figure
    row_height = 5.0
    fig_width = row_height * aspect
    fig = plt.figure(figsize=(fig_width, row_height * len(TARGET_CROPS)))
    gs = gridspec.GridSpec(len(TARGET_CROPS), 1, figure=fig, hspace=0.05)

    # 3. Load Mask
    lr_path = os.path.join(DATA_DIR, "2_StudyArea", basin, "low_runoff_mask.nc")
    ds_lr = xr.open_dataset(lr_path)
    
    last_im = None

    for i, crop_name in enumerate(TARGET_CROPS):
        ax = fig.add_subplot(gs[i], projection=ccrs.PlateCarree())
        crop_files = CROP_MAP[crop_name]
        
        # Calculate Change
        b = get_aggregated_production(BASE_DIR, basin, crop_files)
        s = get_aggregated_production(SCEN_DIR, basin, crop_files)

        # Always plot boundary
        gdf_boundary.boundary.plot(ax=ax, color='black', linewidth=0.8, zorder=3)

        if b is not None and s is not None:
            change = (100 * (s - b) / b.where(b > 0))
            
            # Step 1: Changes in Production
            im = change.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=prod_cmap, norm=prod_norm, add_colorbar=False, zorder=2)
            
            # Step 2: Low Flow Grey Mask
            lr_mask = ds_lr["Low_Runoff"].reindex_like(change, method='nearest').where(change.notnull())
            lr_mask.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=grey_cmap, add_colorbar=False, zorder=1,alpha=0.5)
            
            last_im = im

        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.axis('off')

        # Scale bar on the bottom plot
        if i == len(TARGET_CROPS) - 4:
            width_deg = lon_max - lon_min
            s_len = 500 if width_deg > 10 else 200
            add_scale_bar(ax, s_len, basin)

    # 4. Shared Colorbar
    if last_im:
        cbar_ax = fig.add_axes([0.2, 0.07, 0.6, 0.015])
        cbar = fig.colorbar(last_im, cax=cbar_ax, orientation='horizontal')
        cbar.set_ticks(prod_bounds)
        cbar.ax.set_xticklabels([str(b) for b in prod_bounds], fontsize=20)
        cbar.outline.set_visible(False)

    plt.savefig(os.path.join(FIG_OUT, f"Production_Change_{basin}.png"), dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    for basin in STUDY_AREAS:
        print(f"Processing Basin: {basin}")
        plot_basin_production(basin)