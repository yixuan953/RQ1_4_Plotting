# This code is used to plot the N, P critical losses [kg/ha] through runoff  

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from shapely.ops import unary_union
from matplotlib.gridspec import GridSpec

# === Input and output directories ===
DataDir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method2"
shp_base = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_shp_StudyArea"
PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Boundary_Test/Method2"

Basins = ["Rhine", "Indus", "Yangtze", "LaPlata"]

# Years for statistics
# start_year, end_year = 2005, 2015
year = 2015

# === Loss types and plotting scale ===
loss_types = ["N_runoff", "P_runoff"]
vmins = [0, 0]
vmaxs = [10000000, 1000000]
cmaps = ["YlGnBu", "YlGnBu"]

data_dict = {}

# === Load data ===
for basin in Basins:
    print(f"=== Processing basin: {basin} ===")

    nc_path_N_runoff = f"{DataDir}/{basin}_Total_N_critical_runoff_1986-2015.nc"
    nc_path_P_runoff = f"{DataDir}/{basin}_Total_P_critical_runoff_1986-2015.nc"

    if not (os.path.exists(nc_path_N_runoff) and os.path.exists(nc_path_P_runoff)):
        print(f"Missing data for {basin}, skipping.")
        continue

    # Aggregate to annual critical N, P losses for selected years
    # ds_N = xr.open_dataset(nc_path_N_runoff).sel(Year=slice(f"{start_year}, f"{end_year}"))["N_critical_runoff"]
    # ds_P = xr.open_dataset(nc_path_P_runoff).sel(Year=slice(f"{start_year}",f"{end_year}"))["P_critical_runoff"]
    # N_runoff = ds_N.resample(time="YE").sum(dim="Year", skipna=True)
    # P_runoff = ds_P.resample(time="YE").sum(dim="Year", skipna=True)

    # Annual input 
    ds_N = xr.open_dataset(nc_path_N_runoff).sel(Year=slice(f"{year}"))["N_critical_runoff"]
    ds_P = xr.open_dataset(nc_path_P_runoff).sel(Year=slice(f"{year}"))["P_critical_runoff"]  
    N_runoff = ds_N.where(ds_N != 0)
    P_runoff = ds_P.where(ds_P != 0) 

    avg_N_runoff = N_runoff.mean(dim="Year", skipna=True)
    avg_P_runoff = P_runoff.mean(dim="Year", skipna=True)

    data_dict[("N_runoff", basin)] = avg_N_runoff
    data_dict[("P_runoff", basin)] = avg_P_runoff 

# === Compute basin extents ===
basin_bounds = {}
width_ratios = []
for basin in Basins:
    shp_path = f"{shp_base}/{basin}/{basin}.shp"
    if os.path.exists(shp_path):
        reader = shapereader.Reader(shp_path)
        geoms = list(reader.geometries())
        union_geom = unary_union(geoms)
        xmin, ymin, xmax, ymax = union_geom.bounds
        basin_bounds[basin] = (xmin, ymin, xmax, ymax)
    else:
        xmin, ymin, xmax, ymax = (-180, -90, 180, 90)
        basin_bounds[basin] = (xmin, ymin, xmax, ymax)

    # Width ratio proportional to basin aspect (width / height)
    width = xmax - xmin
    height = ymax - ymin
    width_ratios.append(width / height)

# Normalize width ratios so GridSpec works correctly
total_width = sum(width_ratios)
width_ratios = [w / total_width for w in width_ratios]

# === Plot ===
fig = plt.figure(figsize=(20, 9))
gs = GridSpec(nrows=2, ncols=len(Basins), figure=fig,
              height_ratios=[1, 1], width_ratios=width_ratios)

for i, loss in enumerate(loss_types):
    for j, basin in enumerate(Basins):
        ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())

        if (loss, basin) not in data_dict:
            ax.set_visible(False)
            continue

        da = data_dict[(loss, basin)]
        if "Year" in da.dims and len(da["Year"]) == 1:
            da = da.isel(Year=0)
        im = da.plot.imshow(ax=ax, transform=ccrs.PlateCarree(),
                            cmap=cmaps[i], vmin=vmins[i], vmax=vmaxs[i],
                            add_colorbar=False)

        # Basin boundary
        shp_file = f"{shp_base}/{basin}/{basin}.shp"
        if os.path.exists(shp_file):
            ax.add_geometries(shapereader.Reader(shp_file).geometries(),
                              ccrs.PlateCarree(), facecolor="none",
                              edgecolor="black", linewidth=0.8)

        # River overlay
        river_file = f"{shp_base}/{basin}/{basin}_River.shp"
        if os.path.exists(river_file):
            ax.add_geometries(shapereader.Reader(river_file).geometries(),
                              ccrs.PlateCarree(), facecolor="none",
                              edgecolor="deepskyblue", linewidth=1.0, alpha=0.8)

        ax.set_frame_on(False)
        # --- Set map extent to preserve original aspect ratio ---
        xmin, ymin, xmax, ymax = basin_bounds[basin]
        ax.set_extent([xmin, xmax, ymin, ymax], crs=ccrs.PlateCarree())

        # Titles
        if i == 0:
            ax.set_title(basin, fontsize=18)
        if j == 0:
            ax.text(-0.15, 0.5, loss.replace("_", " "),
                    va="center", ha="right", fontsize=18,
                    rotation=90, transform=ax.transAxes)

        # --- Remove axis ticks ---
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.set_xlabel(''); ax.set_ylabel('')

# === Colorbars ===
cbar_height = 0.02
row_positions = [0.08, 0.52]  # bottom-to-top position for N and P runoff

label_map = {
    "N_runoff": "Critical N runoff",
    "P_runoff": "Critical P runoff"}

for i, loss in enumerate(loss_types):
    cbar_ax = fig.add_axes([0.12, row_positions[1 - i], 0.78, cbar_height])
    sm = plt.cm.ScalarMappable(cmap=cmaps[i])
    sm.set_clim(vmins[i], vmaxs[i])
    cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")

    display_label = label_map.get(loss, loss)  # default to original if not mapped
    cb.set_label(f"{display_label} (kg)", fontsize=15)
    cb.ax.tick_params(labelsize=15)

fig.subplots_adjust(left=0.06, right=0.96, bottom=0.08, top=0.94,
                    wspace=0.03, hspace=0.35)

# === Save ===
out_path = os.path.join(PlotDir, "Critical_NP_Runoff_maincrop_kg.png")
plt.savefig(out_path, dpi=300)
plt.close()

print(f"✅ Figure saved: {out_path}")