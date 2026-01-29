import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import geopandas as gpd
from matplotlib.ticker import ScalarFormatter
import numpy as np

Studyarea =  ["Yangtze"] # ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes =  ["mainrice"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/SM/Fig1"

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
        Basin_mask = ds_input["Basin_mask"].where(ds_input["Total_HA"] > 2500)
        lon = ds_input["lon"].values
        lat = ds_input["lat"].values

        # Variables for crop production plot           
        Rainfed_HA = ds_input["Rainfed_HA"].where(Basin_mask == 1)             
        Irrigated_HA = ds_input["Irrigated_HA"].where(Basin_mask == 1)  
        # Set HA< 6000 to 0 (less than 5% of the pixel)
        Rainfed_HA = Rainfed_HA.where(Rainfed_HA >= 6000, 0)
        Irrigated_HA = Irrigated_HA.where(Irrigated_HA >= 6000, 0)
        total_HA = Rainfed_HA + Irrigated_HA

        # Variables for irrigation plot
        Total_irrigation_amount = ds_input["Total_irrigation_amount"].where(Irrigated_HA > 0)             
        Sus_irrigation_amount = ds_input["Sus_irrigation_amount"].where(Irrigated_HA > 0)             
        Irri_exceedance = ds_input["Irri_exceedance"].where(Irrigated_HA > 0)    
        # Variables for N runoff plot   
        N_Runoff = ds_input["N_Runoff"].where(total_HA > 0)             
        Crit_N_Runoff = ds_input["Crit_N_Runoff"].where(Basin_mask == 1)               
        N_exceedance = ds_input["N_exceedance"].where(total_HA > 0) 
        # Variables for P runoff plot        
        P_Runoff = ds_input["P_Runoff"].where(total_HA > 0)             
        Crit_P_Runoff = ds_input["Crit_P_Runoff"].where(Basin_mask == 1)               
        P_exceedance = ds_input["P_exceedance"].where(total_HA > 0)


        Avg_Yield_Rainfed = ds_input["Avg_Yield_Rainfed"].where(Basin_mask == 1)            
        Avg_Yield_Irrigated = ds_input["Avg_Yield_Irrigated"].where(Basin_mask == 1)
        Total_Production = Rainfed_HA.fillna(0) * Avg_Yield_Rainfed.fillna(0) + Irrigated_HA.fillna(0) * Avg_Yield_Irrigated.fillna(0)
        Total_Production = Total_Production.where(Basin_mask == 1) 

        # ========================== Start the plotting! =========================
        # --------- Figure 1: Total irrigaiton v.s. sustainble irrigation [m3]
        vmin = 0
        vmax = 1000000000
        fig1 = plt.figure(figsize=(14, 5)) 
        gs = GridSpec(1, 3, figure=fig1, width_ratios=[1, 1, 0.6])
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
        ax2.set_title(f"Critical Irrigation Amount")
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
        max_height = 0.7 
        if sus_irrigation_frac / 100 < 1:
            baseline_height = max_height
            sus_height = max_height * sus_irrigation_frac / 100
        else:
            sus_height = max_height
            baseline_height = max_height * 100 / sus_irrigation_frac
        # Plotting the two bars
        ax3.set_ylim(0, 1.0)
        ax3.set_yticks([])
        ax3.bar(0, baseline_height, 0.4, label="Baseline", color="#1376bdf3")
        ax3.bar(1, sus_height, 0.4, label="Critical",  color="#bcdef7f7")
        ax3.set_xlim(-0.6, 1.6)
        # Set X-axis ticks and labels
        ax3.set_xticks([0, 1])
        ax3.set_xticklabels(["Baseline", "Critical"], fontsize=10)
        # Clean up the spines (optional: removes top and right lines for a modern look)
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        # Text for Baseline (Left Bar)
        baseline_str = f"{basin_total_irrigation / 10**9:.2f} × $10^{{9}}$"
        ax3.text(0, baseline_height + 0.025, f"{baseline_str} $m^{{3}}$", ha="center", va="bottom", fontsize=9)
        # Text for Critical (Right Bar)
        sus_total_amount = basin_total_irrigation * (sus_irrigation_frac / 100)
        sus_total_str = f"{sus_total_amount / 10**9:.2f} × $10^{{9}}$"
        ax3.text(1, sus_height + 0.025, f"{sus_total_str} $m^{{3}}$", ha="center", va="bottom", fontsize=9)

        fig1.savefig(os.path.join(output_dir,f"{basin}_{crop}_Irri.png"))
        print(f"Successfully plotted {basin}_{crop}_Irri.png")

        # --------- Figure 2: Total N Runoff v.s. Critical N Runoff [ktons]
        vmin = 0
        vmax = 5
        fig2 = plt.figure(figsize=(12, 5)) # Indus
        gs = GridSpec(1, 3, figure=fig2, width_ratios=[1, 1, 0.6])
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
        ax1.set_title(f"Baseline N losses through cropland runoff")
        ax1.axis('off')        
        # Map 2
        ax2 = fig2.add_subplot(gs[0, 1], projection=proj)
        Crit_N_Runoff = Crit_N_Runoff*0.000001 # Transform to ktons
        im2 = Crit_N_Runoff.plot(ax=ax2,cmap="Reds",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax2.set_title(f"Critical N losses through cropland runoff")
        ax2.axis('off')
        # Shared color bar
        cbar_ax = fig2.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
        cbar = fig2.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("ktons")    

        # Bar Plot
        # --- N Runoff Bar Plot ---
        ax3 = fig2.add_subplot(gs[0, 2]) # Assuming this is ax4 on fig2
        max_height_n = 0.7 
        n_frac = basin_total_crit_N_loss_frac

        if n_frac / 100 < 1:
            n_baseline_height = max_height_n
            n_crit_height = max_height_n * n_frac / 100
        else:
            n_crit_height = max_height_n
            n_baseline_height = max_height_n * 100 / n_frac

        # Plotting the two bars
        ax3.bar(0, n_baseline_height, 0.4, label="Baseline", color="#d65327")
        ax3.bar(1, n_crit_height, 0.4, label="Critical", color="#e5a08a")

        # Configure Axes
        ax3.set_ylim(0, 1.0)
        ax3.set_xlim(-0.6, 1.6)
        ax3.set_yticks([]) # Remove Y ticks as requested
        ax3.set_xticks([0, 1])
        ax3.set_xticklabels(["Baseline", "Critical"], fontsize=10)

        # Clean up spines
        for spine in ['top', 'right', 'left']:
            ax3.spines[spine].set_visible(False)

        # Text for Baseline N (Left Bar)
        # basin_total_N_loss / 10**6 for ktons (assuming input is kg)
        n_baseline_val = basin_total_N_loss / 10**6
        ax3.text(0, n_baseline_height + 0.025, f"{n_baseline_val:.2f} ktons", 
                ha="center", va="bottom", fontsize=9)

        # Text for Critical N (Right Bar)
        n_crit_total = n_baseline_val * (n_frac / 100)
        ax3.text(1, n_crit_height + 0.025, f"{n_crit_total:.2f} ktons", 
                ha="center", va="bottom", fontsize=9)

        fig2.savefig(os.path.join(output_dir,f"{basin}_{crop}_N_runoff.png"))
        print(f"Successfully plotted {basin}_{crop}_N_runoff.png")

        # --------- Figure 3: Total P Runoff v.s. Critical P Runoff [ktons]
        vmin = 0
        vmax = 0.2
        fig3 = plt.figure(figsize=(12, 5)) # Indus
        gs = GridSpec(1, 3, figure=fig3, width_ratios=[1, 1, 0.6])
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
        ax1.set_title(f"Baseline P losses through cropland runoff")
        ax1.axis('off')        
        # Map 2
        ax2 = fig3.add_subplot(gs[0, 1], projection=proj)
        Crit_P_Runoff = Crit_P_Runoff*0.000001 # Transform to ktons
        im2 = Crit_P_Runoff.plot(ax=ax2,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
        basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
        river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
        ax2.set_title(f"Critical P losses through cropland runoff")
        ax2.axis('off')
        # Shared color bar
        cbar_ax = fig3.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
        cbar = fig3.colorbar(im1, cax=cbar_ax, orientation="horizontal")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
        cbar.ax.xaxis.set_major_formatter(formatter)
        cbar.set_label("ktons")           

        # Bar Plot
        # --- N Runoff Bar Plot ---
        ax3 = fig3.add_subplot(gs[0, 2]) # Assuming this is ax4 on fig2
        max_height_P = 0.7 
        P_frac = basin_total_crit_P_loss_frac

        if P_frac / 100 < 1:
            P_baseline_height = max_height_P
            P_crit_height = max_height_P * P_frac / 100
        else:
            P_crit_height = max_height_P
            P_baseline_height = max_height_P * 100 / P_frac

        # Plotting the two bars
        ax3.bar(0, P_baseline_height, 0.4, label="Baseline",  color = "#fbc104" )
        ax3.bar(1, P_crit_height, 0.4, label="Critical", color="#f4db91")

        # Configure Axes
        ax3.set_ylim(0, 1.0)
        ax3.set_xlim(-0.6, 1.6)
        ax3.set_yticks([]) # Remove Y ticks as requested
        ax3.set_xticks([0, 1])
        ax3.set_xticklabels(["Baseline", "Critical"], fontsize=10)

        # Clean up spines
        for spine in ['top', 'right', 'left']:
            ax3.spines[spine].set_visible(False)

        # Text for Baseline P (Left Bar)
        P_baseline_val = basin_total_P_loss / 10**6
        ax3.text(0, P_baseline_height + 0.025, f"{P_baseline_val:.2f} ktons", 
                ha="center", va="bottom", fontsize=9)

        # Text for Critical P (Right Bar)
        P_crit_total = P_baseline_val * (P_frac / 100)
        ax3.text(1, P_crit_height + 0.025, f"{P_crit_total:.2f} ktons", 
                ha="center", va="bottom", fontsize=9)

        fig3.savefig(os.path.join(output_dir,f"{basin}_{crop}_P_runoff.png"))
        print(f"Successfully plotted {basin}_{crop}_P_runoff.png")