# This script is used to plot how much N, P fertilizer input has been changed in each basin, and where the changes are spatially located

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
Studyarea = ["LaPlata", "Rhine", "Indus", "Yangtze"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/1_Fertilizer_input"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig2/2a"
os.makedirs(output_dir, exist_ok=True)

Variables = {
    "N": {"label": "N fertilizer input change (kg N/ha)", "cmap": "PiYG", "vmin": -150, "vmax": 150},
    "P": {"label": "P fertilizer input change (kg P/ha)", "cmap": "PiYG", "vmin": -30, "vmax": 30}
}

def plot_spatial_change(var_key, crop_name, data_dict, basins, output_path):
    # 1. Calculate spans for consistent map heights across basins
    max_lat_span = 0
    basin_info = {}
    for basin in basins:
        gdf = gpd.read_file(os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp"))
        minx, miny, maxx, maxy = gdf.total_bounds
        basin_info[basin] = {'bounds': (minx, miny, maxx, maxy), 
                             'lat_span': maxy - miny, 'lon_span': maxx - minx}
        if (maxy - miny) > max_lat_span: max_lat_span = maxy - miny

    # 2. Layout Setup: Width ratios based on longitude span
    ratios = [basin_info[b]['lon_span'] for b in basins]
    fig = plt.figure(figsize=(26, 8)) 
    gs = gridspec.GridSpec(1, 4, width_ratios=ratios, figure=fig, wspace=0.5)
    
    # 3. Colormap and Normalization (Diverging for changes)
    v_meta = Variables[var_key]
    norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=v_meta["vmin"], vmax=v_meta["vmax"])
    cmap = plt.get_cmap(v_meta["cmap"])

    for i, basin in enumerate(basins):
        ax = fig.add_subplot(gs[i], projection=ccrs.PlateCarree())
        ax.spines['geo'].set_visible(False)
        
        # 4. Background Features
        ax.add_feature(cfeature.LAND, facecolor="#f5f5f5f2")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor="gray")
        
        # 5. Low Runoff Mask (Grey overlay with zorder)
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
        if os.path.exists(low_runoff_path):
            with xr.open_dataset(low_runoff_path) as ds_lr:
                lr_mask = ds_lr["Low_Runoff"] 
            lr_mask.plot.imshow(ax=ax, transform=ccrs.PlateCarree(), 
                                cmap=mcolors.ListedColormap(["#C4C4C4FF"]), 
                                add_colorbar=False, zorder=2)

        # 6. Plot Main Data
        if basin in data_dict:
            data = data_dict[basin]
            im = data.plot.imshow(
                ax=ax, transform=ccrs.PlateCarree(),
                norm=norm, cmap=cmap,
                add_colorbar=False, add_labels=False,zorder=4
            )
            ax.set_title("")
        
        # 7. Basin and River Boundary
        basin_shp = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp")
        
        if os.path.exists(basin_shp):
            gpd.read_file(basin_shp).boundary.plot(ax=ax, color='black', linewidth=1.2, zorder=5)

        # 8. Set View Extent
        minx, miny, maxx, maxy = basin_info[basin]['bounds']
        center_y = (maxy + miny) / 2
        view_miny, view_maxy = center_y - (max_lat_span/2) - 1, center_y + (max_lat_span/2) + 1
        ax.set_extent([minx - 1, maxx + 1, view_miny, view_maxy], crs=ccrs.PlateCarree())
        
        # 9. Ticks/Coordinates
        gl = ax.gridlines(draw_labels=True, linewidth=0, color='none')
        gl.top_labels = False; gl.right_labels = False
        gl.xlocator = MultipleLocator(5); gl.ylocator = MultipleLocator(5)
        gl.xlabel_style = {'size': 14}; gl.ylabel_style = {'size': 14}

    # 10. Vertical Colorbar on the Right
    cbar_ax = fig.add_axes([0.93, 0.25, 0.012, 0.5]) 
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', extend='both')
    cbar.outline.set_visible(False)
    cbar.set_label(v_meta["label"], fontsize=18, labelpad=15)
    cbar.ax.tick_params(labelsize=14)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()

# --- Main Processing Loop ---
for crop in Croptypes:
    data_dict_N = {}
    data_dict_P = {}
    
    for basin in Studyarea:
        input_file = os.path.join(input_dir, f"{basin}_{crop}_Fert_Change.nc")
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
        
        if not (os.path.exists(input_file) and os.path.exists(low_runoff_path)):
            continue
            
        with xr.open_dataset(low_runoff_path) as ds_lr:
            mask_valid = xr.where(ds_lr["Low_Runoff"].isnull(), 1, np.nan)
        
        with xr.open_dataset(input_file) as ds_in:
            data_dict_N[basin] = ds_in["N_change_kgperha"] * mask_valid
            data_dict_P[basin] = ds_in["P_change_kgperha"] * mask_valid

    # Plot N for this crop across 4 basins
    if data_dict_N:
        out_N = os.path.join(output_dir, f"AllBasins_{crop}_N_Change.png")
        plot_spatial_change("N", crop, data_dict_N, Studyarea, out_N)
        
    # Plot P for this crop across 4 basins
    if data_dict_P:
        out_P = os.path.join(output_dir, f"AllBasins_{crop}_P_Change.png")
        plot_spatial_change("P", crop, data_dict_P, Studyarea, out_P)

    print(f"Finished multi-basin maps for {crop}")