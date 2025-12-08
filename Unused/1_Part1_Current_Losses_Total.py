# This code is used to calculate the current losses for main crops and plot the excessive losses [kg]

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

# Output direction
Irrigated_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Output_Unsus_Irrigation"
Rainfed_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Output_Rainfed"

Basins = ["Rhine", "Indus", "Yangtze", "LaPlata"]
Crops = ["winterwheat", "mainrice", "secondrice", "maize", "soybean"]
# start_year, end_year = 2005, 2015
year = 1985

# === Loss types and plotting scale ===
loss_types = ["N_runoff", "P_runoff"]
vmins = [-10000000, -1000000]
vmaxs = [10000000, 1000000]
cmaps = ["PuOr", "PuOr"]

data_dict = {}

# === Load data ===
for basin in Basins:
    print(f"=== Processing basin: {basin} ===")

    # Initialize totals for this basin
    total_N_loss = None
    total_P_loss = None
    total_HA = None

    for crop in Crops:
        irrigated_nc = os.path.join(Irrigated_dir , f"{basin}_{crop}_annual.nc")
        if not os.path.exists(irrigated_nc):
            print(f"Missing file: {irrigated_nc}, skipping crop.")
            continue
        ds_irrigated = xr.open_dataset(irrigated_nc).sel(time=year)
        actual_N_runoff_irrigated = ds_irrigated["N_surf"] + ds_irrigated["N_sub"]
        actual_P_runoff_irrigated = ds_irrigated["P_surf"] + ds_irrigated["P_sub"]

        rainfed_nc = os.path.join(Rainfed_dir , f"{basin}_{crop}_annual.nc")
        if not os.path.exists(rainfed_nc):
            print(f"Missing file: {rainfed_nc}, skipping crop.")
            continue
        ds_rainfed = xr.open_dataset(rainfed_nc).sel(time=year)
        actual_N_runoff_rainfed = ds_rainfed["N_surf"] + ds_rainfed["N_sub"]
        actual_P_runoff_rainfed = ds_rainfed["P_surf"] + ds_rainfed["P_sub"]

        # Harvested area
        HA_rainfed = xr.open_dataset(os.path.join(DataDir, basin, "Mask", f"{crop}_HA_rainfed.nc"))
        HA_rainfed = HA_rainfed["area_total"] if "area_total" in HA_rainfed else list(HA_rainfed.data_vars.values())[0]

        HA_irrigated = xr.open_dataset(os.path.join(DataDir, basin, "Mask", f"{crop}_HA_Irrigated.nc"))
        HA_irrigated = HA_irrigated["area_total"] if "area_total" in HA_irrigated else list(HA_irrigated.data_vars.values())[0]

        # Ensure matched grid
        HA_rainfed = HA_rainfed.interp_like(actual_N_runoff_irrigated)
        HA_irrigated = HA_irrigated.interp_like(actual_N_runoff_irrigated)

        # Compute total HA and loss for this crop
        total_HA_crop = (HA_rainfed.fillna(0) + HA_irrigated.fillna(0))
        total_N_loss_crop = actual_N_runoff_rainfed.fillna(0) * HA_rainfed.fillna(0) + \
                            actual_N_runoff_irrigated.fillna(0) * HA_irrigated.fillna(0)
        total_P_loss_crop = actual_P_runoff_rainfed.fillna(0) * HA_rainfed.fillna(0) + \
                            actual_P_runoff_irrigated.fillna(0) * HA_irrigated.fillna(0)

        # Accumulate across crops
        if total_N_loss is None:
            total_N_loss = total_N_loss_crop
            total_P_loss = total_P_loss_crop
            total_HA = total_HA_crop
        else:
            total_N_loss = total_N_loss + total_N_loss_crop
            total_P_loss = total_P_loss + total_P_loss_crop
            total_HA = total_HA + total_HA_crop

    # Load the total allowable runoff
    nc_path_N_runoff = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method1/{basin}_crit_N_runoff_kg.nc"
    nc_path_P_runoff = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method1/{basin}_crit_P_runoff_kg.nc"

    if not (os.path.exists(nc_path_N_runoff) and os.path.exists(nc_path_P_runoff)):
        print(f"Missing data for {basin}, skipping.")
        continue
    # Annual input 
    ds_N = xr.open_dataset(nc_path_N_runoff).sel(year=slice(f"{year}"))["critical_total_maincrop_N_runoff"]
    ds_P = xr.open_dataset(nc_path_P_runoff).sel(year=slice(f"{year}"))["critical_total_maincrop_P_runoff"]  
    N_runoff = ds_N
    P_runoff = ds_P  

    mask = xr.where(total_HA > 2500, 1, np.nan)

    # Calculate the excess and save it to .nc file
    avg_N_runoff = mask * N_runoff.mean(dim="year", skipna=True)
    avg_P_runoff = mask * P_runoff.mean(dim="year", skipna=True)

    exceed_N_runoff = total_N_loss - avg_N_runoff
    exceed_P_runoff = total_P_loss - avg_P_runoff

    exceed_pixels = xr.where((exceed_N_runoff >0) & (exceed_P_runoff >0), 1 , 0)
    exceed_pixels = exceed_pixels.rename("exceed_pixels")

    # Optionally add metadata
    exceed_pixels.attrs["long_name"] = "Pixels where N and P runoff exceed average"
    exceed_pixels.attrs["units"] = "-"  # dimensionless mask
    exceed_pixels.to_netcdf(f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method1/{basin}_exceed_pixels.nc")
    print("Saved: exceed_pixels.nc")

    data_dict[("N_runoff", basin)] = exceed_N_runoff
    data_dict[("P_runoff", basin)] = exceed_P_runoff 

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

label_map = {
    "N_runoff": "Excessive N runoff",
    "P_runoff": "Excessive P runoff"}

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
out_path = os.path.join(PlotDir, "Excessive_NP_Runoff_maincrop_kg.png")
plt.savefig(out_path, dpi=300)
plt.close()

print(f"✅ Figure saved: {out_path}")