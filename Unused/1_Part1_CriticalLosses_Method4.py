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
PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/0_Test/Boundary_Test/Method4"
BaseDataPath = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/Test_CriticalNP/Method4"

Basins = ["Rhine", "Indus", "Yangtze", "LaPlata"]
crop_map = {"Rice": "RICE", "Wheat": "WHEA", "Maize": "MAIZ", "Soybean": "SOYB"}
loss_types = ["N_runoff", "P_runoff"]

# Define plotting specs (Matching your template format)
unit_config = {
    "kgperha": {
        "suffix": "kgperha.nc",
        "label": "kg/ha",
        "vmins": [0, 0],
        "vmaxs": [25, 2.5],
        "cmaps": ["BuPu", "BuPu"]
    },
    "ktons": {
        "suffix": "kg.nc",
        "label": "ktons",
        "vmins": [0, 0],
        "vmaxs": [0.5, 0.02],
        "cmaps": ["YlGnBu", "YlGnBu"]
    }
}

# === Start Processing ===
for crop, ha_prefix in crop_map.items():
    for unit_key, config in unit_config.items():
        data_dict = {}
        basin_bounds = {}
        width_ratios = []
        active_basins_in_loop = []

        # === Load data ===
        for basin in Basins:
            ha_file = f"{DataDir}/{basin}/Harvest_Area/{ha_prefix}_Harvest_Area_05d_{basin}.nc"
            range_file = f"{DataDir}/{basin}/range.nc"
            nc_N = f"{BaseDataPath}/{crop}/{basin}_crit_N_runoff_{config['suffix']}"
            nc_P = f"{BaseDataPath}/{crop}/{basin}_crit_P_runoff_{config['suffix']}"

            if not (os.path.exists(ha_file) and os.path.exists(nc_N) and os.path.exists(nc_P)):
                continue

            # Load Data
            ds_N = xr.open_dataset(nc_N)
            ds_P = xr.open_dataset(nc_P)
            
            # Setup Masking
            range_mask = xr.open_dataset(range_file)["mask"]
            ha_mask_raw = xr.open_dataset(ha_file)["Harvest_Area"]
            ha_mask = xr.where(ha_mask_raw > 0, 1, np.nan)
            mask = range_mask * ha_mask
            
            N_val = ds_N["Boundary_N_runoff"]
            P_val = ds_P["Boundary_P_runoff"]

            # Unit Conversion (kg to ktons) if necessary
            if unit_key == "ktons":
                N_val *= 1e-6
                P_val *= 1e-6

            data_dict[("N_runoff", basin)] = N_val * mask
            data_dict[("P_runoff", basin)] = P_val * mask
            active_basins_in_loop.append(basin)

            # === Compute basin extents ===
            shp_path = f"{shp_base}/{basin}/{basin}.shp"
            if os.path.exists(shp_path):
                reader = shapereader.Reader(shp_path)
                union_geom = unary_union(list(reader.geometries()))
                xmin, ymin, xmax, ymax = union_geom.bounds
                basin_bounds[basin] = (xmin, ymin, xmax, ymax)
            else:
                xmin, ymin, xmax, ymax = (-180, -90, 180, 90)
                basin_bounds[basin] = (xmin, ymin, xmax, ymax)

            width = xmax - xmin
            height = ymax - ymin
            width_ratios.append(width / height)

        if not active_basins_in_loop:
            continue

        # Normalize width ratios for GridSpec
        total_width_ratio = sum(width_ratios)
        norm_ratios = [w / total_width_ratio for w in width_ratios]

        # === Plot ===
        fig = plt.figure(figsize=(len(active_basins_in_loop) * 5, 9))
        gs = GridSpec(nrows=2, ncols=len(active_basins_in_loop), figure=fig,
                      height_ratios=[1, 1], width_ratios=norm_ratios)

        for i, loss in enumerate(loss_types):
            for j, basin in enumerate(active_basins_in_loop):
                ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())

                if (loss, basin) not in data_dict:
                    ax.set_visible(False)
                    continue

                da = data_dict[(loss, basin)]
                im = da.plot.imshow(ax=ax, transform=ccrs.PlateCarree(),
                                    cmap=config["cmaps"][i], 
                                    vmin=config["vmins"][i], 
                                    vmax=config["vmaxs"][i],
                                    add_colorbar=False)

                # Basin boundary
                shp_file = f"{shp_base}/{basin}/{basin}.shp"
                if os.path.exists(shp_file):
                    ax.add_geometries(shapereader.Reader(shp_file).geometries(),
                                      ccrs.PlateCarree(), facecolor="none",
                                      edgecolor="black", linewidth=0.8)

                # River overlay (Optional, included based on your template)
                river_file = f"{shp_base}/{basin}/{basin}_River.shp"
                if os.path.exists(river_file):
                    ax.add_geometries(shapereader.Reader(river_file).geometries(),
                                      ccrs.PlateCarree(), facecolor="none",
                                      edgecolor="deepskyblue", linewidth=1.0, alpha=0.8)

                ax.set_frame_on(False)
                xmin, ymin, xmax, ymax = basin_bounds[basin]
                ax.set_extent([xmin, xmax, ymin, ymax], crs=ccrs.PlateCarree())

                # Titles and Labels
                if i == 0:
                    ax.set_title(basin, fontsize=18)
                if j == 0:
                    ax.text(-0.15, 0.5, loss.replace("_", " "),
                            va="center", ha="right", fontsize=18,
                            rotation=90, transform=ax.transAxes)

                ax.set_xticks([]); ax.set_yticks([])
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.set_xlabel(''); ax.set_ylabel('')

        # === Colorbars (Exactly as per your template) ===
        cbar_height = 0.02
        row_positions = [0.08, 0.52]  # bottom-to-top position

        for i, loss in enumerate(loss_types):
            cbar_ax = fig.add_axes([0.12, row_positions[1 - i], 0.78, cbar_height])
            sm = plt.cm.ScalarMappable(cmap=config["cmaps"][i])
            sm.set_clim(config["vmins"][i], config["vmaxs"][i])
            cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
            cb.set_label(f"{loss.replace('_', ' ')} ({config['label']})", fontsize=14)
            cb.ax.tick_params(labelsize=12)

        fig.subplots_adjust(left=0.06, right=0.96, bottom=0.08, top=0.94,
                            wspace=0.03, hspace=0.35)

        # === Save ===
        os.makedirs(PlotDir, exist_ok=True)
        out_path = os.path.join(PlotDir, f"Boundary_Loss_{crop}_{unit_key}.png")
        plt.savefig(out_path, dpi=300)
        plt.close()

        print(f"✅ Figure saved: {out_path}")