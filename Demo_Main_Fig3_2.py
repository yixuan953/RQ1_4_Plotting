# This code is used to plot N, P exceedance (avg. 2010 - 2020)

import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import geopandas as gpd
from matplotlib.ticker import ScalarFormatter
import numpy as np

Studyarea =  ["Indus", "LaPlata", "Yangtze", "Rhine"]

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs"

for basin in Studyarea:
    basin_shp_file = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp")
    basin_gdf = gpd.read_file(basin_shp_file)
    river_shp_file = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}_River.shp")
    river_gdf = gpd.read_file(river_shp_file)

    input_file = os.path.join(input_dir, f"{basin}_all_crop_summary.nc") 
    if not os.path.exists(input_file):
        print(f"{basin} basin does not have summary.nc")
        continue

    ds_input = xr.open_dataset(input_file)
    lon = ds_input["lon"].values
    lat = ds_input["lat"].values

    # Variables for N runoff plot   
    N_Runoff = ds_input["N_Runoff"]            
    Crit_N_Runoff = ds_input["Crit_N_Runoff"]               
    N_exceedance = (N_Runoff - Crit_N_Runoff)/N_Runoff
    # Variables for P runoff plot        
    P_Runoff = ds_input["P_Runoff"]           
    Crit_P_Runoff = ds_input["Crit_P_Runoff"]             
    P_exceedance = (P_Runoff - Crit_P_Runoff)/P_Runoff

    # ========== Start plotting ============
    vmin = -100
    vmax = 100
    fig1 = plt.figure(figsize=(10, 6)) 
    gs = GridSpec(1, 2, figure=fig1, width_ratios=[1, 1]) 
    proj = ccrs.PlateCarree()   

    # Map 1
    ax1 = fig1.add_subplot(gs[0, 0], projection=proj)
    N_exceedance = N_exceedance*100
    im1 = N_exceedance.plot(ax=ax1,cmap="PRGn_r",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax1.set_title(f"N exceedance/N critical")
    ax1.axis('off')        
    # Map 2
    ax2 = fig1.add_subplot(gs[0, 1], projection=proj)
    P_exceedance = P_exceedance*100
    im2 = P_exceedance.plot(ax=ax2,cmap="PRGn_r",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax2.set_title(f"P exceedance/P critical")
    ax2.axis('off')

    # Shared color bar
    cbar_ax = fig1.add_axes([0.20, 0.12, 0.60, 0.04])  # [left, bottom, width, height]
    cbar = fig1.colorbar(im1, cax=cbar_ax, orientation="horizontal")
    formatter = ScalarFormatter(useMathText=True)
    # formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
    cbar.ax.xaxis.set_major_formatter(formatter)
    cbar.set_label("%")    

    fig1.savefig(os.path.join(output_dir,f"{basin}_NP_exceedance.png"))
    print(f"Successfully plotted {basin}_NP_exceedance.png")    