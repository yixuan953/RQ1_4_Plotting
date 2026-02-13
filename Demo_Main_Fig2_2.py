# This code is used to plot how much N, P input has been changed

import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import geopandas as gpd
import numpy as np

# Config
Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean", "totalrice"]

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/1_Fertilizer_input"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs/Fig2"

# Standard RdYlGn: Negative=Red, Positive=Green
cmap_change = plt.cm.RdYlGn

for basin in Studyarea:
    basin_shp = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp")
    river_shp = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}_River.shp")
    
    if not os.path.exists(basin_shp): continue
    basin_gdf = gpd.read_file(basin_shp).to_crs("EPSG:4326")
    river_gdf = gpd.read_file(river_shp).to_crs("EPSG:4326")

    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)

    for crop in Croptypes:
        input_file = os.path.join(input_dir, f"{basin}_{crop}_Fert_Change.nc") 
        if not os.path.exists(input_file): continue  

        ds_input = xr.open_dataset(input_file) 
        N_change = ds_input["N_change_kgperha"] * mask_not_low_runoff
        P_change = ds_input["P_change_kgperha"] * mask_not_low_runoff 

        fig = plt.figure(figsize=(14, 7)) 
        gs = GridSpec(1, 2, figure=fig)
        proj = ccrs.PlateCarree()

        # --- Subplot 1: N changes ---
        ax1 = fig.add_subplot(gs[0, 0], projection=proj)
        im1 = N_change.plot(ax=ax1, cmap=cmap_change, transform=proj, 
                            vmin=-150, vmax=150, add_colorbar=True,
                            cbar_kwargs={'label': 'kg N/ha', 'orientation': 'horizontal', 'pad': 0.08})
        ax1.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.plot(ax=ax1, edgecolor="black", linewidth=1, facecolor="none", zorder=10)
        river_gdf.plot(ax=ax1, edgecolor="royalblue", linewidth=0.5, facecolor="none", zorder=11)
        ax1.set_title(f"Nitrogen Change")
        ax1.axis('off')

        # --- Subplot 2: P changes ---
        ax2 = fig.add_subplot(gs[0, 1], projection=proj)
        im2 = P_change.plot(ax=ax2, cmap=cmap_change, transform=proj, 
                            vmin=-30, vmax=30, add_colorbar=True,
                            cbar_kwargs={'label': 'kg P/ha', 'orientation': 'horizontal', 'pad': 0.08})
        ax2.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.plot(ax=ax2, edgecolor="black", linewidth=1, facecolor="none", zorder=10)
        river_gdf.plot(ax=ax2, edgecolor="royalblue", linewidth=0.5, facecolor="none", zorder=11)
        ax2.set_title(f"Phosphorus Change")
        ax2.axis('off')

        plt.suptitle(f"{basin} - {crop.capitalize()}", fontsize=14)
        
        save_path = os.path.join(output_dir, f"{basin}_{crop}_Simple_Map.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Finished: {basin} {crop}")