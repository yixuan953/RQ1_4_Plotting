import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import geopandas as gpd
from matplotlib.ticker import ScalarFormatter
import numpy as np

Studyarea =  ["Indus"] # ["Indus", "LaPlata", "Yangtze", "Rhine"]

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
# input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs"

for basin in Studyarea:
    basin_shp_file = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp")
    basin_gdf = gpd.read_file(basin_shp_file)
    river_shp_file = os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}_River.shp")
    river_gdf = gpd.read_file(river_shp_file)

    input_file = os.path.join(input_dir, f"{basin}_all_crop_summary_baseline.nc") 
    if not os.path.exists(input_file):
        print(f"{basin} basin does not have summary.nc")
        continue

    ds_input = xr.open_dataset(input_file)
    lon = ds_input["lon"].values
    lat = ds_input["lat"].values

    # Mask the low flow 
    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, f"low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_low_runoff:
        low_runoff = ds_low_runoff["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)

    # Variables for N runoff plot   
    N_Runoff = ds_input["N_Runoff"] * mask_not_low_runoff         
    Crit_N_Runoff = ds_input["Crit_N_Runoff"] * mask_not_low_runoff                
    N_exceedance = N_Runoff - Crit_N_Runoff
    # Variables for P runoff plot        
    P_Runoff = ds_input["P_Runoff"] * mask_not_low_runoff            
    Crit_P_Runoff = ds_input["Crit_P_Runoff"] * mask_not_low_runoff              
    P_exceedance =  P_Runoff - Crit_P_Runoff
    # Variables for irrigation      
    Irrigation = ds_input["Irrigation"]           
    Crit_Irrigation = ds_input["Crit_Irrigation"]             
    Irrigation_exceedance =  Irrigation - Crit_Irrigation

    # --------- Figure 1: Total N Runoff v.s. Critical N Runoff [ktons]
    vmin = 0
    vmax = 2
    fig1 = plt.figure(figsize=(12, 5)) # Indus
    gs = GridSpec(1, 3, figure=fig1, width_ratios=[1, 1, 0.6])
    proj = ccrs.PlateCarree()            
    basin_total_N_loss = float(N_Runoff.sum(skipna=True).values)
    basin_total_crit_N_loss = float(Crit_N_Runoff.sum(skipna=True).values)
    # Calculate % exceedance
    if basin_total_N_loss != 0:
        basin_total_crit_N_loss_frac = 100 * basin_total_crit_N_loss / basin_total_N_loss
    else:
        basin_total_crit_N_loss_frac = 0
    # Map 1
    ax1 = fig1.add_subplot(gs[0, 0], projection=proj)
    N_Runoff = N_Runoff*0.000001 # Transform to ktons
    im1 = N_Runoff.plot(ax=ax1,cmap="Reds",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
    ax1.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax1.set_title(f"N losses through cropland runoff")
    ax1.axis('off')        
    # Map 2
    ax2 = fig1.add_subplot(gs[0, 1], projection=proj)
    Crit_N_Runoff = Crit_N_Runoff*0.000001 # Transform to ktons
    im2 = Crit_N_Runoff.plot(ax=ax2,cmap="Reds",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
    ax2.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax2.set_title(f"Critical N losses through cropland runoff")
    ax2.axis('off')
    # Shared color bar
    cbar_ax = fig1.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
    cbar = fig1.colorbar(im1, cax=cbar_ax, orientation="horizontal")
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
    cbar.ax.xaxis.set_major_formatter(formatter)
    cbar.set_label("ktons")    

    # Bar Plot
    # --- N Runoff Bar Plot ---
    ax3 = fig1.add_subplot(gs[0, 2]) # Assuming this is ax4 on fig2
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

    fig1.savefig(os.path.join(output_dir,f"{basin}_total_N_runoff.png"))
    print(f"Successfully plotted {basin}_total_N_runoff.png")

    # --------- Figure 2: Total P Runoff v.s. Critical P Runoff [ktons]
    vmin = 0
    vmax = 0.2
    fig2 = plt.figure(figsize=(12, 5)) # Indus
    gs = GridSpec(1, 3, figure=fig2, width_ratios=[1, 1, 0.6])
    proj = ccrs.PlateCarree()            
    basin_total_P_loss = float(P_Runoff.sum(skipna=True).values)
    basin_total_crit_P_loss = float(Crit_P_Runoff.sum(skipna=True).values)
    # Calculate % exceedance
    if basin_total_P_loss != 0:
         basin_total_crit_P_loss_frac = 100 * basin_total_crit_P_loss / basin_total_P_loss
    else:
        basin_total_crit_P_loss_frac = 0
    # Map 1
    ax1 = fig2.add_subplot(gs[0, 0], projection=proj)
    P_Runoff = P_Runoff*0.000001 # Transform to ktons
    im1 = P_Runoff.plot(ax=ax1,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
    ax1.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax1.set_title(f"P losses through cropland runoff")
    ax1.axis('off')  

    # Map 2
    ax2 = fig2.add_subplot(gs[0, 1], projection=proj)
    Crit_P_Runoff = Crit_P_Runoff*0.000001 # Transform to ktons
    im2 = Crit_P_Runoff.plot(ax=ax2,cmap="YlOrBr",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
    ax2.contourf(low_runoff.lon, low_runoff.lat, low_runoff, colors="#f6efd0da", hatches=['///'], transform=ccrs.PlateCarree())
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax2.set_title(f"Critical P losses through cropland runoff")
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

    fig2.savefig(os.path.join(output_dir,f"{basin}_total_P_runoff.png"))
    print(f"Successfully plotted {basin}_total_P_runoff.png")


    # --------- Figure 3: Irrigation v.s. Critical Irrigation [m3] --> 
    vmin = 0
    vmax = 1
    fig3 = plt.figure(figsize=(12, 5)) # Indus
    gs = GridSpec(1, 3, figure=fig3, width_ratios=[1, 1, 0.6])
    proj = ccrs.PlateCarree()            
    basin_total_Irrigation = float(Irrigation.sum(skipna=True).values)
    basin_total_crit_Irrigation = float(Crit_Irrigation.sum(skipna=True).values)
    # Calculate % exceedance
    if basin_total_Irrigation != 0:
         basin_total_crit_Irrigation_frac = 100 * basin_total_crit_Irrigation / basin_total_Irrigation
    else:
        basin_total_crit_Irrigation_frac = 0

    # Map 1
    ax1 = fig3.add_subplot(gs[0, 0], projection=proj)
    Irrigation = Irrigation / 10**9
    im1 = Irrigation.plot(ax=ax1,cmap="Blues",transform=ccrs.PlateCarree(),add_colorbar=False, vmin = vmin, vmax=vmax)
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax1,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax1.set_title(f"Irrigation Amount")
    ax1.axis('off')  

    # Map 2
    ax2 = fig3.add_subplot(gs[0, 1], projection=proj)
    Crit_Irrigation = Crit_Irrigation / 10**9
    im2 = Crit_Irrigation.plot(ax=ax2,cmap="Blues",transform=ccrs.PlateCarree(),add_colorbar=False,vmin = vmin, vmax=vmax)
    basin_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="black",linewidth=1.5,facecolor="none",zorder=10)
    river_gdf.to_crs("EPSG:4326").plot(ax=ax2,edgecolor="blue",linewidth=0.8,facecolor="none",zorder=11)
    ax2.set_title(f"Critical Irrigation Amount")
    ax2.axis('off')
    # Shared color bar
    cbar_ax = fig3.add_axes([0.20, 0.12, 0.40, 0.04])  # [left, bottom, width, height]
    cbar = fig3.colorbar(im1, cax=cbar_ax, orientation="horizontal")
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 2))  # Show sci when <1e-2 or >1e2
    cbar.ax.xaxis.set_major_formatter(formatter)
    cbar.set_label(f"$10^{{9}}$ m3")           

    # Bar Plot
    # --- Irrigation Bar Plot ---
    ax3 = fig3.add_subplot(gs[0, 2]) # Assuming this is ax4 on fig2
    max_height_Irrigation = 0.7 
    Irrigation_frac = basin_total_crit_Irrigation_frac

    if Irrigation_frac / 100 < 1:
        Irrigation_baseline_height = max_height_Irrigation
        Irrigation_crit_height = max_height_Irrigation * Irrigation_frac / 100
    else:
        Irrigation_crit_height = max_height_Irrigation
        Irrigation_baseline_height = max_height_Irrigation * 100 / Irrigation_frac

    # Plotting the two bars
    ax3.bar(0, Irrigation_baseline_height, 0.4, label="Baseline",  color = "#1f77b4" )
    ax3.bar(1, Irrigation_crit_height, 0.4, label="Critical", color = "#c9e3f6")

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
    Irrigation_baseline_val = basin_total_Irrigation / 10**9
    ax3.text(0, Irrigation_baseline_height + 0.025, f"{Irrigation_baseline_val:.2f} × $10^{{9}}$ m3", 
            ha="center", va="bottom", fontsize=9)

    # Text for Critical P (Right Bar)
    Irrigation_crit_total = Irrigation_baseline_val * (Irrigation_frac / 100)
    ax3.text(1, Irrigation_crit_height + 0.025, f"{Irrigation_crit_total:.2f} × $10^{{9}}$ m3", 
            ha="center", va="bottom", fontsize=9)

    fig3.savefig(os.path.join(output_dir,f"{basin}_total_irri.png"))
    print(f"Successfully plotted {basin}_total_irri.png")
