import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
from matplotlib.colors import ListedColormap, BoundaryNorm


# --- 1. Configuration ---
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
elements = ["N", "P"]

input_dir = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "4_RQ1_Analysis_Results/V5_Plots/SM_Fig/"
    "Share_Agri_Runoff/ncFiles"
)

data_dir = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "2_RQ1_Data"
)

fig_base_dir = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "4_RQ1_Analysis_Results/V5_Plots/SM_Fig/"
    "Share_Agri_Runoff"
)

os.makedirs(fig_base_dir, exist_ok=True)


# --- 2. Color Schemes & Config ---
grey_cmap = ListedColormap(["#EBE9E9"])

n_colors = [
    "#F8B2B2",
    "#E4A5AD",
    "#CF98A9",
    "#BB8BA4",
    "#A67E9F",
    "#92719B",
    "#7D6496",
    "#695791",
    "#544A8D",
    "#403D88",
]

n_boundaries = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
n_cmap = ListedColormap(n_colors)
n_norm = BoundaryNorm(n_boundaries, n_cmap.N)
n_major_ticks = [0, 20, 40, 60, 80, 100]


p_colors = [
    "#F8B2B2",
    "#E4A5AD",
    "#CF98A9",
    "#BB8BA4",
    "#A67E9F",
    "#92719B",
    "#7D6496",
    "#695791",
    "#544A8D",
    "#403D88",
]

p_boundaries = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
p_cmap = ListedColormap(p_colors)
p_norm = BoundaryNorm(p_boundaries, p_cmap.N)
p_major_ticks = [0, 20, 40, 60, 80, 100]


# Use N and P as keys because plot_total_basin receives var_type="N" or "P".
VAR_CONFIG = {
    "N": {
        "cmap": n_cmap,
        "norm": n_norm,
        "ticks": n_major_ticks,
        "tick_labels": [str(t) for t in n_major_ticks],
    },
    "P": {
        "cmap": p_cmap,
        "norm": p_norm,
        "ticks": p_major_ticks,
        "tick_labels": [str(t) for t in p_major_ticks],
    },
}


# --- 3. Helper Functions ---
def add_scale_bar(ax, length_km, basin):
    """
    Add a scale bar to a Cartopy map.

    The map extent must be set before this function is called.
    """
    lon0, lon1, lat0, lat1 = ax.get_extent(
        crs=ccrs.PlateCarree()
    )

    center_lat = (lat0 + lat1) / 2

    geod = cgeo.Geodesic()

    distance_result = geod.inverse(
        np.array([[lon0, center_lat]]),
        np.array([[lon0 + 1, center_lat]]),
    )

    dist_1deg = distance_result[0, 0]

    if not np.isfinite(dist_1deg) or dist_1deg <= 0:
        print(f"Warning: could not calculate scale bar for {basin}.")
        return

    bar_width_deg = (length_km * 1000) / dist_1deg

    if basin == "Yangtze":
        x_start_frac = 0.02

    elif basin == "Rhine":
        bar_frac = bar_width_deg / (lon1 - lon0)
        x_start_frac = 1.25 - bar_frac

    else:
        bar_frac = bar_width_deg / (lon1 - lon0)
        x_start_frac = 0.95 - bar_frac

    x_start_lon = lon0 + (lon1 - lon0) * x_start_frac
    y_pos_lat = lat0 + (lat1 - lat0) * 0.05

    ax.plot(
        [x_start_lon, x_start_lon + bar_width_deg],
        [y_pos_lat, y_pos_lat],
        transform=ccrs.PlateCarree(),
        color="black",
        linewidth=2.0,
        zorder=10,
        clip_on=False,
    )

    text_x = x_start_lon + bar_width_deg / 2

    ax.text(
        text_x,
        y_pos_lat + (lat1 - lat0) * 0.02,
        f"{length_km} km",
        transform=ccrs.PlateCarree(),
        ha="center",
        va="bottom",
        fontsize=30,
        zorder=10,
    )


def select_first_time_step(data_array):
    """
    Select the first time step only when a time dimension exists.
    """
    if "time" in data_array.dims:
        return data_array.isel(time=0)

    return data_array


def plot_total_basin(basin, var_type):
    """
    Generate one map for one basin and one nutrient element.

    The figure layout remains a single Cartopy map with one horizontal
    colorbar.
    """
    if var_type not in VAR_CONFIG:
        raise ValueError(
            f"Unsupported var_type: {var_type}. "
            f"Expected one of {list(VAR_CONFIG)}."
        )

    color_cfg = VAR_CONFIG[var_type]

    variable_name = f"{var_type}_agri_surface_runoff_share"

    nc_path = os.path.join(
        input_dir,
        f"{basin}_{var_type}_agri_surface_runoff_share_2015.nc",
    )

    shp_path = os.path.join(
        data_dir,
        "2_shp_StudyArea",
        basin,
        f"{basin}.shp",
    )

    low_runoff_path = os.path.join(
        data_dir,
        "2_StudyArea",
        basin,
        "low_runoff_mask.nc",
    )

    # Check required input files before creating the figure.
    missing_files = [
        path
        for path in [nc_path, shp_path, low_runoff_path]
        if not os.path.exists(path)
    ]

    if missing_files:
        print(
            f"Skipping {basin} {var_type}. "
            f"Missing file(s): {missing_files}"
        )
        return

    # Read basin boundary.
    gdf_boundary = gpd.read_file(shp_path)

    if gdf_boundary.empty:
        print(f"Skipping {basin} {var_type}: boundary shapefile is empty.")
        return

    lon_min, lat_min, lon_max, lat_max = gdf_boundary.total_bounds

    lat_range = lat_max - lat_min
    lon_range = lon_max - lon_min

    if lat_range <= 0 or lon_range <= 0:
        print(f"Skipping {basin}: invalid shapefile bounds.")
        return

    aspect = lon_range / lat_range

    # Keep the existing single-panel layout.
    fig_height = 6.0

    fig, ax = plt.subplots(
        figsize=(fig_height * aspect, fig_height),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )

    try:
        # Load the nutrient-share data into memory before closing the file.
        with xr.open_dataset(nc_path) as ds:
            if variable_name not in ds:
                available_variables = list(ds.data_vars)

                raise KeyError(
                    f"Variable '{variable_name}' was not found in "
                    f"{nc_path}. Available variables: "
                    f"{available_variables}"
                )

            share_agri_runoff = select_first_time_step(
                ds[variable_name]
            ).load()

        # Load the low-runoff mask into memory.
        with xr.open_dataset(low_runoff_path) as ds_lr:
            if "Low_Runoff" not in ds_lr:
                available_variables = list(ds_lr.data_vars)

                raise KeyError(
                    f"Variable 'Low_Runoff' was not found in "
                    f"{low_runoff_path}. Available variables: "
                    f"{available_variables}"
                )

            low_runoff = select_first_time_step(
                ds_lr["Low_Runoff"]
            ).load()

        # Set extent before calculating the scale-bar position.
        ax.set_extent(
            [lon_min, lon_max, lat_min, lat_max],
            crs=ccrs.PlateCarree(),
        )

        # Main nutrient-share layer.
        im = share_agri_runoff.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=color_cfg["cmap"],
            norm=color_cfg["norm"],
            add_colorbar=False,
            zorder=1,
        )

        # Match the low-runoff mask to the plotted data grid.
        current_lr_mask = low_runoff.reindex_like(
            share_agri_runoff,
            method="nearest",
        )

        # Show grey only where:
        # 1. the nutrient-share data are valid; and
        # 2. Low_Runoff is greater than zero.
        current_lr_mask = current_lr_mask.where(
            share_agri_runoff.notnull()
            & (current_lr_mask > 0)
        )

        current_lr_mask.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=grey_cmap,
            add_colorbar=False,
            zorder=2,
        )

        # Basin boundary.
        gdf_boundary.boundary.plot(
            ax=ax,
            color="black",
            linewidth=1.0,
            zorder=3,
        )

        # Scale bar.
        bar_len = 500 if lon_range > 10 else 200
        add_scale_bar(ax, bar_len, basin)

        # One colorbar per figure.
        cbar = fig.colorbar(
            im,
            ax=ax,
            orientation="horizontal",
            pad=0.08,
            shrink=0.7,
        )

        cbar.set_ticks(color_cfg["ticks"])

        cbar.ax.set_xticklabels(
            color_cfg["tick_labels"],
            fontsize=20,
        )

        cbar.outline.set_visible(False)

        ax.axis("off")

        out_path = os.path.join(
            fig_base_dir,
            f"{basin}_{var_type}.png",
        )

        fig.savefig(
            out_path,
            dpi=300,
            bbox_inches="tight",
        )

        print(f"Saved: {out_path}")

    except Exception as exc:
        print(f"Error plotting {basin} {var_type}: {exc}")

    finally:
        plt.close(fig)


# --- 4. Run ---
for basin in Studyareas:
    print(f"Processing basin: {basin}")

    for element in elements:
        plot_total_basin(basin, element)

print("Plotting complete.")