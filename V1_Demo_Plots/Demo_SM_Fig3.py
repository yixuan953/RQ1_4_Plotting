# This code is used to show the P runoff and exceedance changes over time:
# Year 2010, 2015, 2020 (The end of 2009, 2014, 2019)

import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import geopandas as gpd
from matplotlib.ticker import ScalarFormatter
import numpy as np

Studyarea =  ["Indus"] # ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes =  ["winterwheat"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

# Input directory 1 - Simulated yield and losses
# Baseline scenario directories
# sce = "Baseline"
# Raifed_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline_rainfed"
# Irrigated_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline"

# Sustainable irrigation scenario directories 
# sce = "Sus_Irri"
# Irrigated_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_2_Sus_Irrigation"  # Sustainable irrigation
# Raifed_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline_rainfed"

# Reduced fertilizer scenario directories
sce = "Red_fert"
Raifed_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed/Red_prop"
Irrigated_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Red_prop"

# Input directory 2 - Model input data, including, Harvested area from SPAM2010, Irrigation from VIC-WUR
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"

# Input directory 3 - Calculated critical N, P runoff
crit_loss_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/2_Critical_NP_losses/Method3"

# Figure output director
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/SM/Fig3"

for basin in Studyarea:

    # Load basin range
    basin_mask_file = os.path.join(data_dir, basin, f"range.nc")
    ds_basin_mask = xr.open_dataset(basin_mask_file)
    mask = ds_basin_mask["mask"]
    # Mask the low flow 
    low_runoff_path = os.path.join(data_dir, basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)

    # shp files for plotting
    basin_shp_file = f"/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_shp_StudyArea/{basin}/{basin}.shp"
    basin_gdf = gpd.read_file(basin_shp_file)
    river_shp_file = f"/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_shp_StudyArea/{basin}/{basin}_River.shp"
    river_gdf = gpd.read_file(river_shp_file)

    for crop in Croptypes:

        # Load harvested area
        if crop == "winterwheat":
            cropname = "WHEA"
        elif crop == "maize":
            cropname = "MAIZ"
        elif crop == "soybean":
            cropname = "SOYB"
        elif crop == "mainrice" and basin != "Yangtze":
            cropname = "RICE"
        elif crop == "mainrice" and basin == "Yangtze":
            cropname = "MAINRICE"
        elif crop == "secondrice":
            cropname = "SECONDRICE" 

        Total_HA_file = os.path.join(data_dir, basin, f"Harvest_Area/{cropname}_Harvest_Area_05d_{basin}.nc")
        Rainfed_HA_file = os.path.join(data_dir, basin, f"Harvest_Area/{cropname}_Rainfed_Harvest_Area_05d_{basin}.nc")
        Irrigated_HA_file = os.path.join(data_dir, basin, f"Harvest_Area/{cropname}_Irrigated_Harvest_Area_05d_{basin}.nc")
        if not os.path.exists (Total_HA_file) or not os.path.exists (Rainfed_HA_file) or not os.path.exists (Irrigated_HA_file):
            print (f"Missing {crop} for basin {basin}")
            continue

        ds_total_HA = xr.open_dataset(Total_HA_file)
        Total_HA = ds_total_HA["Harvest_Area"]

        ds_Rainfed_HA = xr.open_dataset(Rainfed_HA_file)
        Rainfed_HA = ds_Rainfed_HA["Harvest_Area"] 

        ds_Irrigated_HA = xr.open_dataset(Irrigated_HA_file)
        Irrigated_HA = ds_Irrigated_HA["Harvest_Area"] 

        lat = mask["lat"].values
        lon = mask["lon"].values    
        mask = mask.where(Total_HA>2500)
       
        # Load simulated P runoff [kg/ha]
        Rainfed_baseline_file = os.path.join(Raifed_baseline_dir, f"{basin}_{crop}_annual.nc")
        Irrigated_baseline_file = os.path.join (Irrigated_baseline_dir, f"{basin}_{crop}_annual.nc")
        if not os.path.exists (Rainfed_baseline_file) or not os.path.exists (Irrigated_baseline_file):
            print (f"Missing {crop} for basin {basin}")
            continue

        ds_rainfed_sim = xr.open_dataset(Rainfed_baseline_file)
        P_Runoff_Rainfed =  ds_rainfed_sim["P_Runoff"].assign_coords(lat=mask.lat, lon=mask.lon)    # [kg/ha]
        
        ds_Irrigated_sim = xr.open_dataset(Irrigated_baseline_file)
        P_Runoff_Irrigated = ds_Irrigated_sim["P_Runoff"].assign_coords(lat=mask.lat, lon=mask.lon)  # [kg/ha]

        if P_Runoff_Irrigated.size == 0 or P_Runoff_Rainfed.size == 0:
                    print(f"Skipping {basin} {crop}: Coordinate mismatch resulted in empty data.")
                    continue
        
        # Load critical N, P losses
        if crop == "winterwheat":
            crop_crit_name = "Wheat"
        elif crop == "maize":
            crop_crit_name = "Maize"
        elif crop == "soybean":
            crop_crit_name = "Soybean"
        elif crop == "mainrice":
            crop_crit_name = "Rice"
        elif crop == "secondrice":
            crop_crit_name = "Rice" 

        crit_P_loss_file = os.path.join(crit_loss_dir, f"{crop_crit_name}/{basin}_crit_P_runoff_kg.nc")
        ds_crit_P_loss = xr.open_dataset(crit_P_loss_file)
        Crit_P_Runoff = 0.000001 * ds_crit_P_loss["Critical_P_runoff"]  # [ktons]
        
        # --- Repeat for 2015 ---
        irri_val_2010 = ds_Irrigated_sim["P_Runoff"].sel(year=2009) * Irrigated_HA 
        rain_val_2010 = ds_rainfed_sim["P_Runoff"].sel(year=2009) * Rainfed_HA 
        P_Runoff_total_2010 = (irri_val_2010 + rain_val_2010) * mask * mask_not_low_runoff * 1e-6

        # --- Repeat for 2015 ---
        irri_val_2015 = ds_Irrigated_sim["P_Runoff"].sel(year=2014) * Irrigated_HA 
        rain_val_2015 = ds_rainfed_sim["P_Runoff"].sel(year=2014) * Rainfed_HA 
        P_Runoff_total_2015 = (irri_val_2015 + rain_val_2015)* mask * mask_not_low_runoff * 1e-6

        # --- Repeat for 2020 ---
        irri_val_2020 = ds_Irrigated_sim["P_Runoff"].sel(year=2019) * Irrigated_HA 
        rain_val_2020 = ds_rainfed_sim["P_Runoff"].sel(year=2019) * Rainfed_HA 
        P_Runoff_total_2020 = (irri_val_2020 + rain_val_2020)* mask * mask_not_low_runoff * 1e-6

        # We use .where(Crit_P_Runoff > 0) to avoid dividing by zero or NaN, 
        P_exceedance_2010 = 100 * mask * (P_Runoff_total_2010 - Crit_P_Runoff) / Crit_P_Runoff.where(Crit_P_Runoff > 0)
        P_exceedance_2015 = 100 * mask * (P_Runoff_total_2015 - Crit_P_Runoff) / Crit_P_Runoff.where(Crit_P_Runoff > 0)
        P_exceedance_2020 = 100 * mask * (P_Runoff_total_2020 - Crit_P_Runoff) / Crit_P_Runoff.where(Crit_P_Runoff > 0)              

        # ========== Plotting starts from here! =============
        # Figures on the first row: P runoff changes 
        vmin_P_runoff = 0   # [ktons]
        vmax_P_runoff = 0.1        

        vmin_exceedance = -100 # [%]
        vmax_exceedance = 100

        fig = plt.figure(figsize = (12, 10))
        gs = GridSpec(2, 3, figure=fig)
        proj = ccrs.PlateCarree() 

        # ------------------- P runoff [ktons] -------------------
        # Map 1
        ax1 = fig.add_subplot(gs[0, 0], projection=proj)
        im1 = P_Runoff_total_2010.plot(ax=ax1,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin_P_runoff, vmax=vmax_P_runoff)
        ax1.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax1.set_title(f"2010")
        ax1.axis('off')        
        # Map 2
        ax2 = fig.add_subplot(gs[0, 1], projection=proj)
        im2 = P_Runoff_total_2015.plot(ax=ax2,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin_P_runoff, vmax=vmax_P_runoff)
        ax2.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax2.set_title(f"2015")
        ax2.axis('off')
        # Map 3
        ax3 = fig.add_subplot(gs[0, 2], projection=proj)
        im3 = P_Runoff_total_2020.plot(ax=ax3,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin_P_runoff, vmax=vmax_P_runoff)
        ax3.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax3,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax3,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax3.set_title(f"2020")
        ax3.axis('off')

        # Shared color bar: For P runoff
        cbar_ax = fig.add_axes([0.20, 0.52, 0.65, 0.02])  # [left, bottom, width, height]
        cbar = fig.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("P runoff [ktons]") 

        # ------------------- P exceedance [%] -------------------
        # Map 1
        ax4 = fig.add_subplot(gs[1, 0], projection=proj)
        im4 = P_exceedance_2010.plot(ax=ax4,cmap="PRGn_r",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin_exceedance, vmax=vmax_exceedance)
        ax4.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax4,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax4,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax4.set_title(f"2010")
        ax4.axis('off')        
        # Map 2
        ax5 = fig.add_subplot(gs[1, 1], projection=proj)
        im5 = P_exceedance_2015.plot(ax=ax5,cmap="PRGn_r",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin_exceedance, vmax=vmax_exceedance)
        ax5.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax5,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax5,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax5.set_title(f"2015")
        ax5.axis('off') 
        # Map 3
        ax6 = fig.add_subplot(gs[1, 2], projection=proj)
        im6 = P_exceedance_2020.plot(ax=ax6,cmap="PRGn_r",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin_exceedance, vmax=vmax_exceedance)
        ax6.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax6,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax6,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax6.set_title(f"2020")
        ax6.axis('off') 

        # Shared color bar: For P exceedance
        cbar_ax = fig.add_axes([0.20, 0.12, 0.65, 0.02])  # [left, bottom, width, height]
        cbar = fig.colorbar(im5, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("P runoff exceedance [%]")  

        fig.savefig(os.path.join(output_dir,f"{basin}_{crop}_P_runoff_{sce}.png"))
        print(f"Successfully plotted {basin}_{crop}_P_runoff_{sce}.png")    