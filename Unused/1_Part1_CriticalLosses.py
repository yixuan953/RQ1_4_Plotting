# This code is used to plot the N, P critical losses [kg/ha] through runoff and leaching 

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs


# Input and output directions
DataDir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
shp_base = '/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_shp_StudyArea'
PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Boundary_Test"

Basins = ["Rhine", "Indus",  "Yangtze", "LaPlata"] 

# Years for statistics 
start_year = 1986
end_year = 2015

# === Scale bar used for plotting ===
loss_types = ["N_leaching", "N_runoff", "P_leaching", "P_runoff"]
vmins = [0, 0, 0, 0]
vmaxs = [200, 20, 3, 3]

data_dict = {}

for basin in Basins:
    print(f"=== Processing basin: {basin} ===")

    nc_path_N_leaching = os.path.join(DataDir, f"{basin}/Hydro/{basin}_N_critical_leaching_annual.nc")
    nc_path_N_runoff = os.path.join(DataDir, f"{basin}/Hydro/{basin}_N_critical_runoff_annual.nc")
    nc_path_P_leaching = os.path.join(DataDir, f"{basin}/Hydro/{basin}_P_critical_leaching_annual.nc")
    nc_path_P_runoff = os.path.join(DataDir, f"{basin}/Hydro/{basin}_P_critical_runoff_annual.nc")

    # output_N_leaching = xr.open_dataset(nc_path_N_leaching).sel(time=slice(start_year, end_year))
    # output_N_runoff = xr.open_dataset(nc_path_N_runoff).sel(time=slice(start_year, end_year))
    # output_P_leaching = xr.open_dataset(nc_path_P_leaching).sel(time=slice(start_year, end_year))
    # output_P_runoff = xr.open_dataset(nc_path_P_runoff).sel(time=slice(start_year, end_year))

    output_N_leaching = xr.open_dataset(nc_path_N_leaching).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))
    output_N_runoff   = xr.open_dataset(nc_path_N_runoff).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))
    output_P_leaching = xr.open_dataset(nc_path_P_leaching).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))
    output_P_runoff   = xr.open_dataset(nc_path_P_runoff).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))

    N_leaching = output_N_leaching["OUT_BASEFLOW"]
    N_runoff = output_N_runoff["OUT_RUNOFF"]
    P_leaching = output_P_leaching["OUT_BASEFLOW"]
    P_runoff = output_P_runoff["OUT_RUNOFF"]

    # Load mask for each basin
    HA_path = os.path.join(DataDir, basin, "Mask", f"{basin}_maize_mask.nc")
    if not os.path.exists(HA_path ):
        print(f"Missing mask file: {HA_path}, skipping crop.")
        continue
    HA_nc = xr.open_dataset(HA_path)
    HA_mask = HA_nc["HA"] 
    mask = xr.where(HA_mask > 0, 1, np.nan)

    avg_N_leaching = mask * (N_leaching.mean(dim='time', skipna=True))
    avg_N_runoff = mask * (N_runoff.mean(dim='time', skipna=True))
    avg_P_leaching = mask * (P_leaching.mean(dim='time', skipna=True))
    avg_P_runoff = mask * (P_runoff.mean(dim='time', skipna=True))

    data_dict[("N_leaching", basin)] = avg_N_leaching
    data_dict[("N_runoff", basin)]   = avg_N_runoff
    data_dict[("P_leaching", basin)] = avg_P_leaching
    data_dict[("P_runoff", basin)]   = avg_P_runoff


# === ~~~ Plotting Start ~~~ ===
# --- Compute each basin’s bounding box and aspect-based width ratios ---
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

# --- Create figure and GridSpec ---
fig = plt.figure(figsize=(20, 17))
gs = GridSpec(
    nrows=4, ncols=len(Basins), figure=fig,
    height_ratios=[1,1,1,1],       # equal row height
    width_ratios=width_ratios      # width scaled by basin aspect
)

# --- Loop through losses × basins ---
for i, loss in enumerate(loss_types):
    for j, basin in enumerate(Basins):
        ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())

        # Skip missing data
        if (loss, basin) not in data_dict:
            ax.set_visible(False)
            continue

        # --- Plot data ---
        da = data_dict[(loss, basin)]
        im = da.plot.imshow(
            ax=ax, transform=ccrs.PlateCarree(),
            cmap='BuPu', vmin=vmins[i], vmax=vmaxs[i],
            add_colorbar=False
        )

        # --- Basin boundary ---
        basin_shp_file = f"{shp_base}/{basin}/{basin}.shp"
        if os.path.exists(basin_shp_file):
            ax.add_geometries(
                shapereader.Reader(basin_shp_file).geometries(),
                ccrs.PlateCarree(),
                facecolor='none', edgecolor='black', linewidth=0.8
            )

        # --- River overlay ---
        river_shp_file = f"{shp_base}/{basin}/{basin}_River.shp"
        if os.path.exists(river_shp_file):
            ax.add_geometries(
                shapereader.Reader(river_shp_file).geometries(),
                ccrs.PlateCarree(),
                facecolor='none', edgecolor='deepskyblue',
                linewidth=1.0, alpha=0.8
            )

        # --- Remove subplot border ---
        ax.set_frame_on(False)

        # --- Set map extent to preserve original aspect ratio ---
        xmin, ymin, xmax, ymax = basin_bounds[basin]
        ax.set_extent([xmin, xmax, ymin, ymax], crs=ccrs.PlateCarree())

        # --- Titles and row labels ---
        if i == 0:
            ax.set_title(basin, fontsize=20)
        if j == 0:
            ax.text(-0.15, 0.5, loss.replace("_", " "),
                    va='center', ha='right', fontsize=20,
                    rotation=90, transform=ax.transAxes)

        # --- Remove axis ticks ---
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.set_xlabel(''); ax.set_ylabel('')

# --- Add horizontal colorbars for each row ---
cbar_height = 0.015
row_positions = [0.045, 0.285, 0.525, 0.76]  # bottom-to-top position of each row

for i, loss in enumerate(loss_types):
    cbar_ax = fig.add_axes([0.12, row_positions[3 - i], 0.78, cbar_height])
    sm = plt.cm.ScalarMappable(cmap='BuPu')
    sm.set_clim(vmins[i], vmaxs[i])
    cb = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cb.set_label(f"{loss} (kg/ha)", fontsize=15)
    cb.ax.tick_params(labelsize=15) 

# --- Adjust spacing ---
fig.subplots_adjust(left=0.06, right=0.96, bottom=0.08, top=0.96,
                    wspace=0.03, hspace=0.50)

# --- Save figure ---
out_path = os.path.join(PlotDir, "Critical_NP_Losses_4x4_withRivers_balanced.png")
plt.savefig(out_path, dpi=300)
plt.close()

print(f"✅ Figure saved with equal row height and basin aspect preserved:\n{out_path}")