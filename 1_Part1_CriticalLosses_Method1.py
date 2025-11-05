# Plot N and P critical losses [kg/ha] through runoff: Method 1

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from shapely.ops import unary_union
from matplotlib.gridspec import GridSpec

# === Paths ===
DataDir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
shp_base = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_shp_StudyArea"
PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Boundary_Test/Method1"

Basins = ["Rhine", "Indus", "Yangtze", "LaPlata"]
start_year, end_year = 2005, 2015

# === Loss types and plotting scale ===
loss_types = ["N_runoff", "P_runoff"]
vmins = [0, 0]
vmaxs = [50, 3]
cmaps = ["BuPu", "BuPu"]

data_dict = {}

# === Load data ===
for basin in Basins:
    print(f"=== Processing basin: {basin} ===")

    nc_path_N_runoff = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method1/{basin}_crit_N_runoff.nc"
    nc_path_P_runoff = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method1/{basin}_crit_P_runoff.nc"

    if not (os.path.exists(nc_path_N_runoff) and os.path.exists(nc_path_P_runoff)):
        print(f"Missing data for {basin}, skipping.")
        continue

    # Aggregate to annual critical N, P losses for selected years
    ds_N = xr.open_dataset(nc_path_N_runoff).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-01"))["critical_maincrop_N_runoff"]
    ds_P = xr.open_dataset(nc_path_P_runoff).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-01"))["critical_maincrop_P_runoff"]

    N_runoff = ds_N.resample(time="YE").sum(dim="time", skipna=True)
    P_runoff = ds_P.resample(time="YE").sum(dim="time", skipna=True)

    # Mask
    HA_path = os.path.join(DataDir, basin, "Mask", f"{basin}_maize_mask.nc")
    if not os.path.exists(HA_path):
        print(f"Missing mask file: {HA_path}, skipping crop.")
        continue

    HA_mask = xr.open_dataset(HA_path)["HA"]
    mask = xr.where(HA_mask > 0, 1, np.nan)

    avg_N_runoff = mask * N_runoff.mean(dim="time", skipna=True)
    avg_P_runoff = mask * P_runoff.mean(dim="time", skipna=True)

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

for i, loss in enumerate(loss_types):
    cbar_ax = fig.add_axes([0.12, row_positions[1 - i], 0.78, cbar_height])
    sm = plt.cm.ScalarMappable(cmap=cmaps[i])
    sm.set_clim(vmins[i], vmaxs[i])
    cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cb.set_label(f"{loss} (kg/ha)", fontsize=14)
    cb.ax.tick_params(labelsize=12)

fig.subplots_adjust(left=0.06, right=0.96, bottom=0.08, top=0.94,
                    wspace=0.03, hspace=0.35)

# === Save ===
out_path = os.path.join(PlotDir, "Critical_NP_Runoff_maincrop.png")
plt.savefig(out_path, dpi=300)
plt.close()

print(f"✅ Figure saved: {out_path}")