import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.colors as mcolors
from shapely.geometry import Point

# --- 1. Configuration ---
STUDY_AREAS = ["LaPlata", "Rhine", "Indus", "Yangtze"]
CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}

DATA_DIR = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
BASE_DIR = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
SCEN_DIR = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/4_Inc_fert"
FIG_OUT  = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Demo_Plots/Fig4_Crop_Prod_Red"
os.makedirs(FIG_OUT, exist_ok=True)

# --- 2. Custom Color Setup ---
colors = ["#8b0000", "#fb4343", "#ff8c00", "#f9cb14", "#fae22d", "#f4ff5e", "#66b535"]
bounds = [-100, -80, -50, -30, -20, -10, 0, 100]
custom_cmap = mcolors.ListedColormap(colors)
custom_norm = mcolors.BoundaryNorm(bounds, custom_cmap.N)

# Define a static grey colormap for the low flow area
grey_cmap = mcolors.ListedColormap(["#ece9e9"]) 

# --- 3. Helper Functions ---

def add_scale_bar(ax, length_km, location=(0.70, 0.04)):
    lon0, lon1, lat0, lat1 = ax.get_extent()
    center_lat = (lat0 + lat1) / 2
    deg_per_km = 1 / (111.32 * np.cos(np.deg2rad(center_lat)))
    deg_width = length_km * deg_per_km
    start_lon = lon0 + (lon1 - lon0) * location[0]
    start_lat = lat0 + (lat1 - lat0) * location[1]
    
    ax.plot([start_lon, start_lon + deg_width], [start_lat, start_lat],
            transform=ccrs.PlateCarree(), color='black', linewidth=3, zorder=10)
    
    ax.text(start_lon + deg_width/2, start_lat + (lat1 - lat0) * 0.03, 
            f"{length_km} km", transform=ccrs.PlateCarree(),
            ha='center', va='bottom', fontsize=16, zorder=11)

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

def plot_production_grid():
    # Process data and masks
    results = {basin: {} for basin in STUDY_AREAS}
    low_flow_masks = {} # Dictionary to store grey area masks

    for basin in STUDY_AREAS:
        lr_path = os.path.join(DATA_DIR, "2_StudyArea", basin, "low_runoff_mask.nc")
        with xr.open_dataset(lr_path) as ds_lr:
            # Areas to plot grey: Where Low_Runoff is NOT null (assuming 1 = low flow)
            low_flow_masks[basin] = xr.where(ds_lr["Low_Runoff"].notnull(), 1, np.nan)
            # Areas to plot data: Where Low_Runoff IS null
            valid_runoff_mask = xr.where(ds_lr["Low_Runoff"].isnull(), 1, np.nan)

        for crop, files in CROP_MAP.items():
            b = get_aggregated_production(BASE_DIR, basin, files)
            s = get_aggregated_production(SCEN_DIR, basin, files)
            if b is not None and s is not None:
                # Production change calculation
                change = (100 * (s - b) / b.where(b > 0))
                # Only keep data in high-runoff areas
                results[basin][crop] = change * valid_runoff_mask

    crops = list(CROP_MAP.keys())
    fig, axes = plt.subplots(len(crops), len(STUDY_AREAS), figsize=(24, 20), 
                             subplot_kw={'projection': ccrs.PlateCarree()})

    for c, basin in enumerate(STUDY_AREAS):
        shp_path = os.path.join(DATA_DIR, "2_shp_StudyArea", basin, f"{basin}.shp")
        shp = gpd.read_file(shp_path)
        bbox = shp.total_bounds
        
        width_deg = bbox[2] - bbox[0]
        s_len = 1000 if width_deg > 15 else 500 if width_deg > 5 else 100

        for r, crop in enumerate(crops):
            ax = axes[r, c]
            
            # Step 1: Plot the Grey Low-Flow Mask as a background
            if basin in low_flow_masks:
                low_flow_masks[basin].plot(ax=ax, transform=ccrs.PlateCarree(),
                                           cmap=grey_cmap, add_colorbar=False, zorder=1)
            
            # Step 2: Plot the Crop Change Data
            data = results[basin].get(crop)
            if data is not None:
                im = data.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=custom_cmap, 
                               norm=custom_norm, add_colorbar=False, zorder=2)
                
                # Step 3: Add Scale Bar (First row only)
                # if r == 3:
                    # add_scale_bar(ax, s_len)
            else:
                ax.text(0.5, 0.5, "N/A", transform=ax.transAxes, ha='center', color='grey')
            
            # Step 4: Add Basin Boundary
            shp.boundary.plot(ax=ax, color='black', linewidth=1.2, zorder=3)
            ax.set_extent([bbox[0], bbox[2], bbox[1], bbox[3]], crs=ccrs.PlateCarree())
            
            ax.axis('off')
            if r == 0: ax.set_title(basin, fontsize=24, pad=25, fontweight='bold')
            if c == 0: ax.text(-0.15, 0.5, crop, transform=ax.transAxes, rotation=90, 
                               va='center', fontsize=24, fontweight='bold')

    # Global Colorbar
    cbar_ax = fig.add_axes([0.3, 0.06, 0.4, 0.02])
    cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal', ticks=bounds[1:-1])
    cb.set_label("Crop Production Change [%]", fontsize=20, labelpad=15)
    cb.ax.set_xticklabels([f"{b}" for b in bounds[1:-1]], fontsize=20)

    plt.subplots_adjust(wspace=0.08, hspace=0.08, bottom=0.12)
    plt.savefig(os.path.join(FIG_OUT, "Crop_Grid_NoScaleBar.png"), dpi=300, bbox_inches='tight')
    print("Success! Plot generated with grey low-flow areas.")

if __name__ == "__main__":
    plot_production_grid()