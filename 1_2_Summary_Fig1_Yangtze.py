import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import geopandas as gpd
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import ListedColormap
import numpy as np


Studyarea =  ["Yangtze"] # ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes =  ["winterwheat", "maize", "mainrice", "secondrice", "soybean"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/1_Summary"

for basin in Studyarea:
    basin_shp_file = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp")
    basin_gdf = gpd.read_file(basin_shp_file)
    river_shp_file = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}_River.shp")
    river_gdf = gpd.read_file(river_shp_file)

    for crop in Croptypes:

        input_file = os.path.join(input_dir, f"{basin}_{crop}_summary_baseline.nc") # Basline scenario
        if not os.path.exists(input_file):
            print(f"{basin} basin does not have {crop}")
            continue
        ds_input = xr.open_dataset(input_file)
        Basin_mask = ds_input["Basin_mask"]
        lon = ds_input["lon"].values
        lat = ds_input["lat"].values
        # Variables for irrigation plot
        Total_irrigation_amount = ds_input["Total_irrigation_amount"].where(Basin_mask == 1)             
        Sus_irrigation_amount = ds_input["Sus_irrigation_amount"].where(Basin_mask == 1)             
        Irri_exceedance = ds_input["Irri_exceedance"].where(Basin_mask == 1)    
        # Variables for N runoff plot   
        N_Runoff = ds_input["N_Runoff"].where(Basin_mask == 1)             
        Crit_N_Runoff = ds_input["Crit_N_Runoff"].where(Basin_mask == 1)               
        N_exceedance = ds_input["N_exceedance"].where(Basin_mask == 1) 
        # Variables for P runoff plot        
        P_Runoff = ds_input["P_Runoff"].where(Basin_mask == 1)             
        Crit_P_Runoff = ds_input["Crit_P_Runoff"].where(Basin_mask == 1)               
        P_exceedance = ds_input["P_exceedance"].where(Basin_mask == 1)
        # Variables for crop production plot           
        Rainfed_HA = ds_input["Rainfed_HA"].where(Basin_mask == 1)             
        Irrigated_HA = ds_input["Irrigated_HA"].where(Basin_mask == 1)             
        Avg_Yield_Rainfed = ds_input["Avg_Yield_Rainfed"].where(Basin_mask == 1)            
        Avg_Yield_Irrigated = ds_input["Avg_Yield_Irrigated"].where(Basin_mask == 1)
        Total_Production = Rainfed_HA.fillna(0) * Avg_Yield_Rainfed.fillna(0) + Irrigated_HA.fillna(0) * Avg_Yield_Irrigated.fillna(0)
        Total_Production = Total_Production.where(Basin_mask == 1) 

        # if crop == "secondrice": # Plot an extra summary for secondrice + mainrice
        #     input_file1 = os.path.join(input_dir, f"{basin}_mainrice_baseline_summary.nc")
        #     input_file2 = os.path.join(input_dir, f"{basin}_secondrice_baseline_summary.nc")
        #     if not os.path.exists(input_file1) or not os.path.exists(input_file2):
        #         print(f"{basin} basin does not have {crop}")
        #         continue  
        #     ds_input1 = xr.open_dataset(input_file1)
        #     ds_input2 = xr.open_dataset(input_file2)      

        # else:

        # ========================== Start the plotting! =========================
        # --------- Figure 1: Total irrigaiton v.s. sustainble irrigation [m3]
        vmin = 0
        vmax = 100000000
        fig1 = plt.figure(figsize=(12, 4.5)) 
        gs = GridSpec(1, 3, figure=fig1, width_ratios=[1, 1, 0.3])
        proj = ccrs.PlateCarree()            
        basin_total_irrigation = float(Total_irrigation_amount.sum(skipna=True).values)
        basin_sus_irrigation = float(Sus_irrigation_amount.sum(skipna=True).values)
        # Calculate % sustainable
        if basin_total_irrigation != 0:
            sus_irrigation_frac = 100 * basin_sus_irrigation / basin_total_irrigation
        else:
            sus_irrigation_frac = 0
        # Map 1
        ax1 = fig1.add_subplot(gs[0, 0], projection=proj)
        im1 = Total_irrigation_amount.plot(ax=ax1,cmap="Blues",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax1.set_title(f"Baseline Irrigation Amount")
        ax1.axis('off')        
        # Map 2
        ax2 = fig1.add_subplot(gs[0, 1], projection=proj)
        im2 = Sus_irrigation_amount.plot(ax=ax2,cmap="Blues",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax2.set_title(f"Sustainable Irrigation Amount")
        ax2.axis('off')
        # Shared color bar
        cbar_ax = fig1.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
        cbar = fig1.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("Irrigation Amount (m³)")           
        # Bar Plot
        ax3 = fig1.add_subplot(gs[0, 2])
        sus_height = .8*sus_irrigation_frac/100
        ax3.bar(0, 0.8, 0.1, label="Current (100%)", alpha=0.4, color = "#c9e3f6")
        ax3.bar(0, sus_height, 0.1, label="Sustainable %", alpha=0.7, color = "#1f77b4" )
        ax3.set_ylim(0, 1.0)
        ax3.axis('off')
        # Add absolute amount text on top
        total_str = f"{basin_total_irrigation / 10**10:.2f} × $10^{{10}}$"
        ax3.text(0.0, 0.85,f"Total irrigation amount:\n{total_str} m3",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, 0.70,f"Unsus:\n{100-sus_irrigation_frac:.1f}%",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, sus_height-0.3,f"Sustain:\n{sus_irrigation_frac:.1f}%",ha="center", va="bottom", fontsize=10)

        fig1.savefig(os.path.join(output_dir,f"{basin}_{crop}_Irri.png"))
        print(f"Successfully plotted {basin}_{crop}_Irri.png")

        # --------- Figure 2: Total N Runoff v.s. Critical N Runoff [ktons]
        vmin = 0
        vmax = 10
        fig2 = plt.figure(figsize=(12, 4.5)) # Indus
        gs = GridSpec(1, 3, figure=fig2, width_ratios=[1, 1, 0.3])
        proj = ccrs.PlateCarree()            
        basin_total_N_loss = float(N_Runoff.sum(skipna=True).values)
        basin_total_crit_N_loss = float(Crit_N_Runoff.sum(skipna=True).values)
        # Calculate % exceedance
        if basin_total_N_loss != 0:
            basin_total_crit_N_loss_frac = 100 * basin_total_crit_N_loss / basin_total_N_loss
        else:
            basin_total_crit_N_loss_frac = 0
        # Map 1
        ax1 = fig2.add_subplot(gs[0, 0], projection=proj)
        N_Runoff = N_Runoff*0.000001 # Transform to ktons
        im1 = N_Runoff.plot(ax=ax1,cmap="Reds",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax1.set_title(f"Current N runoff")
        ax1.axis('off')        
        # Map 2
        ax2 = fig2.add_subplot(gs[0, 1], projection=proj)
        Crit_N_Runoff = Crit_N_Runoff*0.000001 # Transform to ktons
        im2 = Crit_N_Runoff.plot(ax=ax2,cmap="Reds",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax2.set_title(f"Critical N runoff")
        ax2.axis('off')
        # Shared color bar
        cbar_ax = fig2.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
        cbar = fig2.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("N Runoff (ktons)")           
        # Bar Plot
        ax3 = fig2.add_subplot(gs[0, 2])
        sus_height = 0.8*basin_total_crit_N_loss_frac/100
        ax3.bar(0, 0.8, 0.1, label="Current (100%)", alpha=0.4, color = "#e5a08a")
        ax3.bar(0, sus_height, 0.1, label="Critical %", alpha=0.7, color = "#d65327" )
        ax3.set_ylim(0, 1.0)
        ax3.axis('off')
        # Add absolute amount text on top
        total_str = f"{basin_total_N_loss/10**6:.2f}" # From kg to ktons
        ax3.text(0.0, 0.85,f"Total N runoff:\n{total_str} ktons",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, 0.55,f"Exceedance:\n{100-basin_total_crit_N_loss_frac:.1f}%",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, sus_height-0.10,f"Critical:\n{basin_total_crit_N_loss_frac:.1f}%",ha="center", va="bottom", fontsize=10)

        fig2.savefig(os.path.join(output_dir,f"{basin}_{crop}_N_runoff.png"))
        print(f"Successfully plotted {basin}_{crop}_N_runoff.png")

        # --------- Figure 3: Total P Runoff v.s. Critical P Runoff [ktons]
        vmin = 0
        vmax = 1
        fig3 = plt.figure(figsize=(12, 4.5)) # Indus
        gs = GridSpec(1, 3, figure=fig3, width_ratios=[1, 1, 0.3])
        proj = ccrs.PlateCarree()            
        basin_total_P_loss = float(P_Runoff.sum(skipna=True).values)
        basin_total_crit_P_loss = float(Crit_P_Runoff.sum(skipna=True).values)
        # Calculate % exceedance
        if basin_total_P_loss != 0:
            basin_total_crit_P_loss_frac = 100 * basin_total_crit_P_loss / basin_total_P_loss
        else:
            basin_total_crit_P_loss_frac = 0
        # Map 1
        ax1 = fig3.add_subplot(gs[0, 0], projection=proj)
        P_Runoff = P_Runoff*0.000001 # Transform to ktons
        im1 = P_Runoff.plot(ax=ax1,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax1.set_title(f"Baseline P runoff")
        ax1.axis('off')        
        # Map 2
        ax2 = fig3.add_subplot(gs[0, 1], projection=proj)
        Crit_P_Runoff = Crit_P_Runoff*0.000001 # Transform to ktons
        im2 = Crit_P_Runoff.plot(ax=ax2,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax2.set_title(f"Critical P runoff")
        ax2.axis('off')
        # Shared color bar
        cbar_ax = fig3.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
        cbar = fig3.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("P Runoff (ktons)")           
        # Bar Plot
        ax3 = fig3.add_subplot(gs[0, 2])
        sus_height = 0.8*basin_total_crit_P_loss_frac/100
        ax3.bar(0, 0.8, 0.1, label="Current (100%)", alpha=0.4, color = "#fde291")
        ax3.bar(0, sus_height, 0.1, label="Critical %", alpha=0.7, color = "#e8ae00" )
        ax3.set_ylim(0, 1.0)
        ax3.axis('off')
        # Add absolute amount text on top
        total_str = f"{basin_total_P_loss/ 10**6:.2f}" # From kg to ktons
        ax3.text(0.0, 0.85,f"Total P runoff:\n{total_str} ktons",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, 0.55,f"Exceedance:\n{100-basin_total_crit_P_loss_frac:.1f}%",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, sus_height-0.10,f"Critical:\n{basin_total_crit_P_loss_frac:.1f}%",ha="center", va="bottom", fontsize=10)

        fig3.savefig(os.path.join(output_dir,f"{basin}_{crop}_P_runoff.png"))
        print(f"Successfully plotted {basin}_{crop}_P_runoff.png")


        # --------- Figure 4: Current crop production + Contribution by unsustainle area 
        vmin = 0
        vmax = 800
        fig4 = plt.figure(figsize=(12, 4.5)) # Indus
        gs = GridSpec(1, 3, figure=fig4, width_ratios=[1, 1, 0.3])
        proj = ccrs.PlateCarree()   
        
        sus_pixels = Basin_mask.copy()
        sus_pixels = xr.where((N_exceedance + P_exceedance + Irri_exceedance) > 3, 0, 1)
        sus_production = Total_Production * sus_pixels
        basin_total_crop_production = float(Total_Production.sum(skipna=True).values)
        basin_sus_crop_prod = float(sus_production.sum(skipna=True).values)

        # Calculate % exceedance
        if basin_total_crop_production != 0:
            basin_sus_crop_prod_frac = 100 * basin_sus_crop_prod / basin_total_crop_production
        else:
            basin_sus_crop_prod_frac = 0
        # Map 1
        ax1 = fig4.add_subplot(gs[0, 0], projection=proj)
        Total_Production = Total_Production*0.000001 # Transform to ktons
        im1 = Total_Production.plot(ax=ax1,cmap="YlGn",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax1.set_title(f"Baseline production")
        ax1.axis('off')    

        # Map 2
        ax2 = fig4.add_subplot(gs[0, 1], projection=proj)
        mask_irri = xr.where(Irri_exceedance == 11, 1, 0)
        mask_n    = xr.where(N_exceedance == 11, 2, 0)
        mask_p    = xr.where(P_exceedance == 11, 4, 0)
        combined_mask = mask_irri + mask_n + mask_p
        combined_mask = combined_mask.where(combined_mask > 0, np.nan)

        colors = [
            "#1f77b4", # 1: Blue (Irri)
            "#d62728", # 2: Red (N)
            "#9467bd", # 3: Purple (Irri + N)
            "#ffbf00", # 4: Yellow (P)
            "#2ca02c", # 5: Green (Irri + P)
            "#ff7f0e", # 6: Orange (N + P)
            "#7f7f7f"  # 7: Grey (All)
        ]
        my_cmap = ListedColormap(colors)

        im = ax2.pcolormesh(lon, lat, combined_mask, transform=ccrs.PlateCarree(),cmap=my_cmap, vmin=0.5, vmax=7.5,alpha=0.8)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax2.set_title(f"Boundary exceedance")
        ax2.axis('off')

        # Shared color bar
        cbar_ax = fig4.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
        cbar = fig4.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-3, 3))  # Show sci when <1e-2 or >1e2
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("Crop production (ktons)")           
        # Bar Plot
        ax3 = fig4.add_subplot(gs[0, 2])
        sus_height = 0.8*basin_sus_crop_prod_frac/100
        ax3.bar(0, 0.8, 0.1, label="Current (100%)", alpha=0.4, color = "#B8B6AE")
        ax3.bar(0, sus_height, 0.1, label="Critical %", alpha=0.7, color = "#83c076" )
        ax3.set_ylim(0, 1.0)
        ax3.axis('off')

        # Add absolute amount text on top
        total_str = f"{basin_total_crop_production/10**10:.2f} × $10^{{4}}$" # Transform kg to ktons
        ax3.text(0.0, 0.85,f"Total crop production:\n{total_str} ktons",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, 0.55,f"Boundary\nexceeded:\n{100-basin_sus_crop_prod_frac:.1f}%",ha="center", va="bottom", fontsize=10)
        ax3.text(0.0, sus_height-0.10,f"Boundary\nmet:\n{basin_sus_crop_prod_frac:.1f}%",ha="center", va="bottom", fontsize=10)

        fig4.savefig(os.path.join(output_dir,f"{basin}_{crop}_prod.png"))
        print(f"Successfully plotted {basin}_{crop}_prod.png")