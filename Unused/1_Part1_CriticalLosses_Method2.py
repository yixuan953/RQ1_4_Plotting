# This code is used to plot the N, P critical losses [kg/ha] through runoff  

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
from matplotlib.gridspec import GridSpec

# === Input and output directories ===
DataDir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method2/Rainfed"
shp_base = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_shp_StudyArea"
MaskDir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Boundary_Test/Method2/Rainfed"

os.makedirs(PlotDir, exist_ok=True)

Basins = ["Rhine", "Indus", "Yangtze", "LaPlata"]
crops = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

# Years for statistics
start_year = 2014
end_year = 2015

# Scale bars for plotting
vmaxs = [50, 2]  # [N runoff, P runoff]

# ===========================================================
for basin in Basins:
    print(f"\n=== Processing basin: {basin} ===")

    basin_shp_file = f"{shp_base}/{basin}/{basin}.shp"
    river_shp_file = f"{shp_base}/{basin}/{basin}_River.shp"

    for crop in crops:
        print(f"  → Crop: {crop}")

        nc_path_N_runoff = os.path.join(DataDir, f"{basin}_{crop}_N_critical_runoff_1986-2015.nc")
        nc_path_P_runoff = os.path.join(DataDir, f"{basin}_{crop}_P_critical_runoff_1986-2015.nc")

        if not (os.path.exists(nc_path_N_runoff) and os.path.exists(nc_path_P_runoff)):
            print(f"     Missing data for {crop}, skipping.")
            continue

        # --- Open and subset data ---
        ds_N = xr.open_dataset(nc_path_N_runoff).sel(Year=slice(start_year, end_year))
        ds_P = xr.open_dataset(nc_path_P_runoff).sel(Year=slice(start_year, end_year))

        # Harmonize dimension names
        ds_N = ds_N.rename({"Lat": "lat", "Lon": "lon"})
        ds_P = ds_P.rename({"Lat": "lat", "Lon": "lon"})

        N_runoff = ds_N["Runoff"]
        P_runoff = ds_P["Runoff"]

        # --- Load and align mask ---
        HA_path = os.path.join(MaskDir, basin, "Mask", f"{basin}_{crop}_mask.nc")
        if not os.path.exists(HA_path):
            print(f"     Missing mask file: {HA_path}, skipping crop.")
            continue

        HA_nc = xr.open_dataset(HA_path)
        HA_mask = HA_nc["HA"]

        # Rename mask dims to match runoff data
        if "Lat" in HA_mask.dims or "Lon" in HA_mask.dims:
            HA_mask = HA_mask.rename({"Lat": "lat", "Lon": "lon"})

        # Interpolate mask if grids differ
        if HA_mask.sizes != N_runoff.isel(Year=0).sizes:
            print(f"     Adjusting mask grid for {basin}-{crop}")
            HA_mask = HA_mask.interp(lat=N_runoff.lat, lon=N_runoff.lon)

        # --- Compute masked means safely ---
        avg_N_runoff = N_runoff.mean(dim="Year", skipna=True).where(HA_mask > 0)
        avg_P_runoff = P_runoff.mean(dim="Year", skipna=True).where(HA_mask > 0)

        # === Plot ===
        fig = plt.figure(figsize=(10, 5))
        gs = GridSpec(1, 2, figure=fig, wspace=0.15)

        titles = ["Critical N runoff (kg N/ha)", "Critical P runoff (kg P/ha)"]
        data_list = [avg_N_runoff, avg_P_runoff]

        for i, (data, title, vmax) in enumerate(zip(data_list, titles, vmaxs)):
            ax = fig.add_subplot(gs[0, i], projection=ccrs.PlateCarree())
            ax.set_title(title, fontsize=13)

            # Extract coordinates and ensure 2D grids
            lat = data["lat"]
            lon = data["lon"]

            if lat.ndim == 1 and lon.ndim == 1:
                lon2d, lat2d = np.meshgrid(lon, lat)
            else:
                lon2d, lat2d = lon.values, lat.values

            print(f"     {basin}-{crop}: data shape={data.shape}, lat={lat.shape}, lon={lon.shape}")

            # Plot
            im = ax.pcolormesh(
                lon2d, lat2d, data.values,
                transform=ccrs.PlateCarree(),
                cmap="BuPu",
                vmin=0, vmax=vmax,
                shading="auto"
            )

            # Add shapefile boundaries
            if os.path.exists(basin_shp_file):
                shp = shapereader.Reader(basin_shp_file)
                for rec in shp.records():
                    geom = rec.geometry
                    ax.add_geometries([geom], ccrs.PlateCarree(),
                                     facecolor="none", edgecolor="black", linewidth=0.8)

            if os.path.exists(river_shp_file):
                river_shp = shapereader.Reader(river_shp_file)
                for rec in river_shp.records():
                    geom = rec.geometry
                    ax.add_geometries([geom], ccrs.PlateCarree(),
                                     facecolor="none", edgecolor="deepskyblue", linewidth=1.0)

            # Add coastlines and borders
            ax.coastlines(resolution="10m", linewidth=0.5)
            ax.add_feature(cfeature.BORDERS, linewidth=0.3)

            # Colorbar
            cbar = plt.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, shrink=0.8)
            cbar.set_label("[kg/ha]", fontsize=10)

        fig.suptitle(f"{basin} – {crop}", fontsize=14, fontweight="bold")
        plt.tight_layout(rect=[0, 0, 1, 0.95])

        # Save
        out_path = os.path.join(PlotDir, f"{basin}_{crop}_Critical_NP_losses.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        print(f"     ✅ Saved: {out_path}")