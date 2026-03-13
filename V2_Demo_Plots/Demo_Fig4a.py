# This script is used to plot how much runoff has been decreased
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator

# --- Configuration ---
# Rice
Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
InputCrops = ["mainrice", "secondrice"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
fig_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4a/Rice"

# Soybean
# Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
# InputCrops = ["soybean"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
# fig_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4a/Soybean"

# # Winterwheat
# Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
# InputCrops = ["winterwheat"]
# fig_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4a/Wheat"

# # Maize
# Studyareas = ["LaPlata", "Rhine", "Indus", "Yangtze"]
# InputCrops = ["maize"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
# fig_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig4/4a/Maize"

os.makedirs(fig_output_dir, exist_ok=True)

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
sc4_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/4_Inc_fert"

Variables = {
    "Crop_Production": {"label": "Yield Reduction (%)", "cmap": "RdYlGn"}
}

def GetCropProd(file_name):
    ds = xr.open_dataset(file_name)
    mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
    Total_HA = ds["Total_HA"] * mask
    Irri_HA = ds["Irrigated_HA"]
    Yield_Irri = ds["Avg_Yield_Irrigated"]
    Rainfed_HA = ds["Rainfed_HA"]
    Yield_Rainfed = ds["Avg_Yield_Rainfed"]
    Total_crop_prod = Irri_HA * Yield_Irri + Rainfed_HA * Yield_Rainfed
    Total_crop_prod = Total_crop_prod * mask
    
    return Total_HA, Total_crop_prod

def Yield_Reduction_Map(var_key, data_dict, basins, output_path):
    # 1. Calculate spans for normalization (consistent height logic)
    max_lat_span = 0
    basin_info = {}
    for basin in basins:
        gdf = gpd.read_file(os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp"))
        minx, miny, maxx, maxy = gdf.total_bounds
        basin_info[basin] = {'bounds': (minx, miny, maxx, maxy), 
                             'lat_span': maxy - miny, 'lon_span': maxx - minx}
        if (maxy - miny) > max_lat_span: max_lat_span = maxy - miny

    # 2. Layout Setup
    num_basins = len(basins) # Dynamic count (will be 1 for Yangtze)
    ratios = [basin_info[b]['lon_span'] for b in basins]
    
    # Adjust figure width dynamically: 7 inches per basin is usually a good rule
    fig = plt.figure(figsize=(7 * num_basins + 2, 10)) 

    # Use num_basins instead of the hard-coded 4
    gs = gridspec.GridSpec(1, num_basins, width_ratios=ratios, figure=fig, wspace=1.0 if num_basins > 1 else 0)
    norm = mcolors.Normalize(vmin=-100, vmax=100)

    for i, basin in enumerate(basins):
        ax = fig.add_subplot(gs[i], projection=ccrs.PlateCarree())
        ax.spines['geo'].set_visible(False)
        
        # 4. Background Features
        ax.add_feature(cfeature.LAND, facecolor="#f5f5f5f2")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor="gray")
        
        # 5. Low Runoff Mask
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
        if os.path.exists(low_runoff_path):
            with xr.open_dataset(low_runoff_path) as ds_lr:
                lr_mask = ds_lr["Low_Runoff"] 
            lr_mask.plot.imshow(ax=ax, transform=ccrs.PlateCarree(), 
                                cmap=mcolors.ListedColormap(["#D3D3D3FF"]), 
                                add_colorbar=False, zorder=2)

        # 6. Plot Main Data
        data = data_dict[basin]
        im = data.plot.imshow(
            ax=ax, transform=ccrs.PlateCarree(),
            norm=norm, cmap="RdYlGn",
            add_colorbar=False, zorder=4
        )
        
        # 7. Basin Boundary
        gdf = gpd.read_file(os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp"))
        gdf.boundary.plot(ax=ax, color='black', linewidth=1.5, zorder=5)

        # 8. Set View Extent
        minx, miny, maxx, maxy = basin_info[basin]['bounds']
        center_y = (maxy + miny) / 2
        view_miny, view_maxy = center_y - (max_lat_span/2) - 2, center_y + (max_lat_span/2) + 2
        ax.set_extent([minx - 2, maxx + 2, view_miny, view_maxy], crs=ccrs.PlateCarree())
        
        # 9. Ticks/Coordinates
        gl = ax.gridlines(draw_labels=True, linewidth=0, color='none')
        gl.top_labels = False; gl.right_labels = False
        gl.xlocator = MultipleLocator(5); gl.ylocator = MultipleLocator(5)
        gl.xlabel_style = {'size': 15, 'color': 'black'}
        gl.ylabel_style = {'size': 15, 'color': 'black'}

    # 10. Vertical Colorbar on the Right
    # [left, bottom, width, height]
    # left=0.92 puts it outside the subplots area
    cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.025]) 
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.outline.set_visible(False)

    # Formatting the colorbar
    cbar.set_label(Variables[var_key]["label"], fontsize=20, labelpad=15)
    cbar.ax.tick_params(labelsize=15)

    # Use bbox_inches='tight' but ensure we don't clip the new colorbar
    plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()

results = { "Crop_Production": {} }

for basin in Studyareas:
    print(f"Processing Basin: {basin}")

    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)

    base_Production_basin = 0.1; sc4_Production_basin= 0; basin_total_HA = 0
    
    found_data = False
    for crop in InputCrops:
        baseline_file_name = os.path.join(baseline_dir, f"{basin}_{crop}_summary.nc")
        sc4_file_name = os.path.join(sc4_dir, f"{basin}_{crop}_summary.nc")

        if not os.path.exists(baseline_file_name) or not os.path.exists(sc4_file_name) :
            continue 
        
        base_HA, base_Production = GetCropProd(baseline_file_name)
        base_Production_basin += base_Production.fillna(0)
        basin_total_HA += base_HA.fillna(0)

        base_HA_nouse, sc4_Production = GetCropProd(sc4_file_name)
        sc4_Production_basin += sc4_Production.fillna(0)
    
    red_Production_runoff = 100 * (sc4_Production_basin - base_Production_basin)/base_Production_basin 
    mask_low_HA = xr.where(basin_total_HA > 2500, 1, np.nan)
    results["Crop_Production"][basin] = red_Production_runoff * mask_not_low_runoff * mask_low_HA

for var_key in Variables:
    out_file = os.path.join(fig_output_dir, f"Reduction_{var_key}.png")
    Yield_Reduction_Map(var_key, results[var_key], Studyareas, out_file)
    print(f"Saved figure for {var_key}")