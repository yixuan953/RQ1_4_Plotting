# This code is to plot:
# 1) Proportion of yield reduction when reducing irrigation and fertilzer [%]
# 2) Total production reduction due to these practices

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader

# Input and output directions
DataDir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
shp_base = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_shp_StudyArea"

Baseline_irrigated_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Output_Unsus_Irrigation"
Baseline_rainfed_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Output_Rainfed"
Sus_irrigated_dir  = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Output_Sus_Irrigation"
Sus_fert_irrigated_dir  = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Output_Sus_Irrigation"
Sus_rainfed_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Output_Rainfed"

Basins = ["Indus", "Rhine", "LaPlata", "Yangtze"] 
CropTypes = ["mainrice", "maize", "winterwheat", "soybean", "secondrice"]

PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Part3"

# Years for statistics 
start_year = 2005
end_year = 2015

for basin in Basins:
    print(f"=== Processing basin: {basin} ===")

    for crop in CropTypes:
        print(f"> {crop}")

        # Load the output of different scenarios
        output_baseline_irrigated_nc = os.path.join(Baseline_irrigated_dir, f"{basin}_{crop}_annual.nc")
        output_Sus_irrigated_nc = os.path.join(Sus_irrigated_dir, f"{basin}_{crop}_annual.nc")      
        output_Sus_fert_irrigated_nc = os.path.join(Sus_fert_irrigated_dir, f"{basin}_{crop}_annual.nc") 
        output_baseline_rainfed_nc = os.path.join(Baseline_rainfed_dir, f"{basin}_{crop}_annual.nc")
        output_Sus_rainfed_nc = os.path.join(Sus_rainfed_dir, f"{basin}_{crop}_annual.nc")
        if not os.path.exists(output_baseline_irrigated_nc) or not os.path.exists(output_Sus_irrigated_nc) or not os.path.exists(output_Sus_fert_irrigated_nc) or not os.path.exists(output_baseline_rainfed_nc) or not os.path.exists(output_Sus_rainfed_nc):
            print(f"Missing {crop} for {basin}")
            continue

        ds_baseline_irrigated = xr.open_dataset(output_baseline_irrigated_nc).sel(time=slice(start_year, end_year))
        ds_Sus_irrigated = xr.open_dataset(output_Sus_irrigated_nc).sel(time=slice(start_year, end_year))
        ds_Sus_fert_irrigated = xr.open_dataset(output_Sus_fert_irrigated_nc).sel(time=slice(start_year, end_year))
        ds_baseline_rainfed = xr.open_dataset(output_baseline_rainfed_nc).sel(time=slice(start_year, end_year))
        ds_Sus_rainfed = xr.open_dataset(output_Sus_rainfed_nc).sel(time=slice(start_year, end_year))

        Baseline_irrigated_yield = ds_baseline_irrigated["Storage"]
        Sus_irrigated_yield = ds_Sus_irrigated["Storage"]
        Sus_fert_irrigated_yield = ds_Sus_fert_irrigated["Storage"]        
        Baseline_rainfed_yield = ds_baseline_rainfed["Storage"]
        Sus_rainfed_yield = ds_Sus_rainfed["Storage"]

        # Load harvested area for rainfed land
        HA_path_rainfed = os.path.join(DataDir, basin, "Mask", f"{crop}_HA_rainfed.nc")
        if not os.path.exists(HA_path_rainfed):
            print(f"Missing mask file: {HA_path_rainfed}, skipping crop.")
            continue
        HA_nc_rainfed = xr.open_dataset(HA_path_rainfed)
        HA_rainfed = HA_nc_rainfed["area_total"] if "area_total" in HA_nc_rainfed else list(HA_nc_rainfed.data_vars.values())[0]

        # Load harvested area for irrigated land
        HA_path_irrigated = os.path.join(DataDir, basin, "Mask", f"{crop}_HA_Irrigated.nc")
        if not os.path.exists(HA_path_irrigated):
            print(f"Missing mask file: {HA_path_irrigated}, skipping crop.")
            continue
        HA_nc_irrigated = xr.open_dataset(HA_path_irrigated)
        HA_irrigated = HA_nc_irrigated["area_total"] if "area_total" in HA_nc_irrigated else list(HA_nc_irrigated.data_vars.values())[0]
        
        # Calculate the total harvest area and get the mask
        HA_total = HA_rainfed.fillna(0) + HA_irrigated.fillna(0)
        HA_total = xr.DataArray(HA_total, dims=("lat", "lon"))
        HA_total = HA_total.interp(lat=Baseline_irrigated_yield.lat, lon=Baseline_irrigated_yield.lon, method="nearest")
        valid_HA = xr.where(HA_total > 2500, 1, 0)

        # Calculate the average annual yield for different scenarios
        Baseline_irrigated_yield_avg = Baseline_irrigated_yield.mean(dim="time", skipna=True) * valid_HA
        Baseline_rainfed_yield_avg = Baseline_rainfed_yield.mean(dim="time", skipna=True) * valid_HA
        Sus_irrigated_yield_avg = Sus_irrigated_yield.mean(dim="time", skipna=True) * valid_HA
        Sus_fert_irrigated_yield_avg = Sus_fert_irrigated_yield.mean(dim="time", skipna=True) * valid_HA
        Sus_rainfed_yield_avg = Sus_rainfed_yield.mean(dim="time", skipna=True) * valid_HA  

        # Calculate the percentage of yield reduction
        Red_irrigated_yield =  100 * (Baseline_irrigated_yield_avg - Sus_fert_irrigated_yield_avg)/Baseline_irrigated_yield_avg
        Red_rainfed_yield =  100 * (Baseline_rainfed_yield_avg - Sus_rainfed_yield_avg)/Baseline_rainfed_yield_avg

        # Calculate the total crop produciton [tons/year] under each scenario
        Baseline_production = 0.001 * (Baseline_irrigated_yield_avg * HA_irrigated + Baseline_rainfed_yield_avg * HA_rainfed).sum(dim=["lat", "lon"], skipna=True)
        Sus_Irri_production = 0.001 * (Sus_irrigated_yield_avg * HA_irrigated + Baseline_rainfed_yield_avg * HA_rainfed).sum(dim=["lat", "lon"], skipna=True)
        Sus_Irri_Fert_production = 0.001 * (Sus_fert_irrigated_yield_avg * HA_irrigated + Sus_rainfed_yield_avg * HA_rainfed).sum(dim=["lat", "lon"], skipna=True)

        # ========== Plotting start ==============
        fig = plt.figure(figsize=(16, 5))
        shp_file = f"{shp_base}/{basin}/{basin}.shp"
        river_file = f"{shp_base}/{basin}/{basin}_River.shp"

        # ---- Shared normalization for both maps ----
        vmin, vmax = 0, 20

        # 1) Map with overlays
        ax1 = fig.add_subplot(1, 3, 1, projection=ccrs.PlateCarree())
        im1 = ax1.pcolormesh(Red_irrigated_yield.lon, Red_irrigated_yield.lat, Red_irrigated_yield,
                            vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), cmap="Reds")
        ax1.set_title("Yield reduction on irrigated land [%]")

        # Add basin boundary
        if os.path.exists(shp_file):
            ax1.add_geometries(shapereader.Reader(shp_file).geometries(),
                            crs=ccrs.PlateCarree(), facecolor="none",
                            edgecolor="black", linewidth=1.2)

        # Add river
        if os.path.exists(river_file):
            ax1.add_geometries(shapereader.Reader(river_file).geometries(),
                            crs=ccrs.PlateCarree(), facecolor="none",
                            edgecolor="deepskyblue", linewidth=1.0, alpha=0.8)


        # 2) Map without overlays, same color scale
        ax2 = fig.add_subplot(1, 3, 2, projection=ccrs.PlateCarree())
        im2 = ax2.pcolormesh(Red_irrigated_yield.lon, Red_irrigated_yield.lat, Red_irrigated_yield,
                            vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), cmap="Reds")
        ax2.set_title("Yield reduction on rainfed land [%]")


        # ---- Shared Colorbar for both maps ----
        cbar = fig.colorbar(im1, ax=[ax1, ax2], orientation="vertical", shrink=0.75, pad=0.03)
        cbar.set_label("Yield reduction [%]")


        # 3) Bar chart of total production
        ax3 = fig.add_subplot(1, 3, 3)
        scenarios = ["Baseline", "Sus-Irrigation", "Sus-Irr+Fert"]
        values = [
            float(Baseline_production),
            float(Sus_Irri_production),
            float(Sus_Irri_Fert_production)
        ]

        bar_colors = ["#f3d277fc", "#92dff5", "#cba4ea"]  
        bars = ax3.bar(scenarios, values, color=bar_colors)
        ax3.set_ylabel("Production (tons/year)")
        ax3.set_title("Total Production Comparison")

        for bar in bars:
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f"{bar.get_height():.1f}", ha="center", va="bottom")

        # Save figure
        os.makedirs(PlotDir, exist_ok=True)
        outname = os.path.join(PlotDir, f"{basin}_{crop}_Yield_Prod_red.png")
        plt.tight_layout()
        plt.savefig(outname, dpi=300)
        plt.close()

        print(f"Saved: {outname}")