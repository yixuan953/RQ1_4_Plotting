import os

import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import cartopy.crs as ccrs
import cartopy.geodesic as cgeo


# ============================================================
# 1. CONFIGURATION
# ============================================================

STUDY_AREAS = [
    "LaPlata",
    "Indus",
    "Yangtze",
    "Rhine",
]

TARGET_CROPS = [
    "Wheat",
    "Maize",
    "Rice",
    "Soybean",
]

# Displayed crop name: crop file names
CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice"],
    "Soybean": ["soybean"],
}


DATA_DIR = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "2_RQ1_Data"
)

BASE_DIR = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "3_RQ1_Model_Outputs/"
    "4_Analysis4Plotting/"
    "0_Summary/"
    "1_Baseline"
)

SCEN_DIR = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "3_RQ1_Model_Outputs/"
    "4_Analysis4Plotting/"
    "0_Summary/"
    "0_Yp"
)

FIG_OUT = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "4_RQ1_Analysis_Results/"
    "V5_Plots/"
    "SM_Fig/"
    "Yield"
)

os.makedirs(FIG_OUT, exist_ok=True)


# Three yield modes × four basins = twelve figures.
YIELD_MODES = {
    "Irrigated": {
        "folder": BASE_DIR,
        "components": [("Avg_Yield_Irrigated", "Irrigated_HA")],
    },
    "Rainfed": {
        "folder": BASE_DIR,
        "components": [("Avg_Yield_Rainfed", "Rainfed_HA")],
    },
    "Potential": {
        "folder": SCEN_DIR,
        "components": [
            ("Avg_Yield_Irrigated", "Irrigated_HA"),
            ("Avg_Yield_Rainfed", "Rainfed_HA"),
        ],
    },
}


# Minimum harvested area required for plotting a grid cell
MIN_TOTAL_HA = 2500


# ============================================================
# 2. COLOR SETTINGS
# ============================================================

# Green sequential palette. The first two classes are deliberately
# medium green so values from 0 to 2,000 remain clearly visible.
yield_colors = [
    "#B7DFA8",
    "#8FCB82",
    "#68B66B",
    "#49A35B",
    "#2F8F4E",
    "#1F7A43",
    "#146637",
    "#0B522D",
    "#063F24",
    "#032D1A",
]

yield_bounds = [
    0,
    1000,
    2000,
    3000,
    4000,
    5000,
    6000,
    7000,
    8000,
    9000,
    10000,
]

yield_cmap = mcolors.ListedColormap(yield_colors)

yield_norm = mcolors.BoundaryNorm(
    yield_bounds,
    yield_cmap.N,
)

grey_cmap = mcolors.ListedColormap(
    ["#ece9e9"]
)


# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

def add_scale_bar(ax, length_km, basin):
    """
    Add a geographical scale bar to a Cartopy axis.

    Parameters
    ----------
    ax : cartopy.mpl.geoaxes.GeoAxes
        Cartopy map axis.
    length_km : float
        Scale-bar length in kilometres.
    basin : str
        Basin name, used to adjust the scale-bar position.
    """

    lon0, lon1, lat0, lat1 = ax.get_extent(
        crs=ccrs.PlateCarree()
    )

    center_lat = (lat0 + lat1) / 2.0

    geod = cgeo.Geodesic()

    # Distance represented by one degree longitude at center latitude
    distance_result = geod.inverse(
        (lon0, center_lat),
        (lon0 + 1.0, center_lat),
    )

    distance_one_degree = distance_result[0, 0]

    if distance_one_degree <= 0:
        return

    bar_width_deg = (
        length_km * 1000.0
    ) / distance_one_degree

    longitude_span = lon1 - lon0
    latitude_span = lat1 - lat0

    bar_fraction = bar_width_deg / longitude_span

    # Basin-specific scale-bar positions
    if basin in ["Yangtze", "Rhine"]:
        # Bottom-left
        x_start_fraction = 0.04
    else:
        # Bottom-right
        x_start_fraction = 0.96 - bar_fraction

    # Keep the scale bar within the map
    x_start_fraction = max(
        0.02,
        min(x_start_fraction, 0.98 - bar_fraction),
    )

    x_start_lon = (
        lon0 + longitude_span * x_start_fraction
    )

    y_position_lat = (
        lat0 + latitude_span * 0.05
    )

    x_end_lon = x_start_lon + bar_width_deg

    # Scale-bar line
    ax.plot(
        [x_start_lon, x_end_lon],
        [y_position_lat, y_position_lat],
        transform=ccrs.PlateCarree(),
        color="black",
        linewidth=2.0,
        zorder=10,
        clip_on=False,
    )

    # Small end marks
    tick_height = latitude_span * 0.008

    ax.plot(
        [x_start_lon, x_start_lon],
        [
            y_position_lat - tick_height,
            y_position_lat + tick_height,
        ],
        transform=ccrs.PlateCarree(),
        color="black",
        linewidth=2.0,
        zorder=10,
        clip_on=False,
    )

    ax.plot(
        [x_end_lon, x_end_lon],
        [
            y_position_lat - tick_height,
            y_position_lat + tick_height,
        ],
        transform=ccrs.PlateCarree(),
        color="black",
        linewidth=2.0,
        zorder=10,
        clip_on=False,
    )

    # Scale-bar text
    text_x = (
        x_start_lon + bar_width_deg / 2.0
    )

    ax.text(
        text_x,
        y_position_lat + latitude_span * 0.018,
        f"{length_km} km",
        transform=ccrs.PlateCarree(),
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=24,
        color="black",
        zorder=11,
    )


def get_aggregated_yield(
    folder,
    basin,
    crop_list,
    components,
):
    """Calculate area-weighted yield for one crop and yield mode."""

    total_production = None
    total_harvested_area = None

    for crop in crop_list:
        input_path = os.path.join(folder, f"{basin}_{crop}_summary.nc")

        if not os.path.exists(input_path):
            print(f"    Warning: file not found:\n    {input_path}")
            continue

        with xr.open_dataset(input_path) as ds:
            required_variables = ["Basin_mask"]
            for yield_variable, area_variable in components:
                required_variables.extend([yield_variable, area_variable])

            missing_variables = [
                variable for variable in required_variables if variable not in ds
            ]
            if missing_variables:
                print(
                    f"    Warning: missing variables in "
                    f"{os.path.basename(input_path)}: {missing_variables}"
                )
                continue

            basin_mask = ds["Basin_mask"]
            crop_production = None
            crop_area = None

            for yield_variable, area_variable in components:
                area = ds[area_variable].fillna(0)
                production = ds[yield_variable].fillna(0) * area
                crop_production = production if crop_production is None else crop_production + production
                crop_area = area if crop_area is None else crop_area + area

            valid_cells = (
                basin_mask.notnull()
                & (basin_mask > 0)
                & crop_area.notnull()
                & (crop_area > MIN_TOTAL_HA)
            )

            crop_production = crop_production.where(valid_cells).load()
            crop_area = crop_area.where(valid_cells).load()

        if total_production is None:
            total_production = crop_production
            total_harvested_area = crop_area
        else:
            total_production = total_production.fillna(0) + crop_production.fillna(0)
            total_harvested_area = (
                total_harvested_area.fillna(0) + crop_area.fillna(0)
            )

    if total_production is None:
        return None

    valid_total_area = total_harvested_area.where(total_harvested_area > 0)
    aggregated_yield = (total_production / valid_total_area).where(
        valid_total_area.notnull()
    )
    aggregated_yield.name = "Yield"
    return aggregated_yield


def load_low_runoff_mask(basin):
    """
    Load the low-runoff mask for one basin.
    """

    low_runoff_path = os.path.join(
        DATA_DIR,
        "2_StudyArea",
        basin,
        "low_runoff_mask.nc",
    )

    if not os.path.exists(low_runoff_path):
        print(
            f"    Warning: low-runoff mask not found:\n"
            f"    {low_runoff_path}"
        )
        return None

    with xr.open_dataset(low_runoff_path) as ds:

        if "Low_Runoff" not in ds:
            print(
                f"    Warning: variable 'Low_Runoff' "
                f"not found in:\n"
                f"    {low_runoff_path}"
            )
            return None

        low_runoff = ds["Low_Runoff"].load()

    return low_runoff


def prepare_low_runoff_mask(
    low_runoff,
    yield_data,
):
    """
    Align the low-runoff mask to the yield grid.

    Only low-runoff cells containing valid yield data are retained.
    """

    if low_runoff is None:
        return None

    try:
        aligned_mask = low_runoff.reindex_like(
            yield_data,
            method="nearest",
        )
    except (ValueError, KeyError):
        # interp_like can work better when coordinate spacing differs
        aligned_mask = low_runoff.interp_like(
            yield_data,
            method="nearest",
        )

    plotted_mask = xr.where(
        (aligned_mask > 0)
        & yield_data.notnull(),
        1.0,
        np.nan,
    )

    return plotted_mask


# ============================================================
# 4. MAIN PLOTTING FUNCTION
# ============================================================

def plot_basin_yield(
    basin,
    yield_mode,
    mode_config,
):
    """
    Create one yield figure for one basin and one yield mode.

    Each figure contains four vertical panels:
        1. Wheat
        2. Maize
        3. Rice
        4. Soybean
    """

    print(
        f"  Plotting {yield_mode.lower()} yield"
    )

    # --------------------------------------------------------
    # Read basin boundary
    # --------------------------------------------------------

    shapefile_path = os.path.join(
        DATA_DIR,
        "2_shp_StudyArea",
        basin,
        f"{basin}.shp",
    )

    if not os.path.exists(shapefile_path):
        print(
            f"    Error: shapefile not found:\n"
            f"    {shapefile_path}"
        )
        return

    boundary_gdf = gpd.read_file(
        shapefile_path
    )

    if boundary_gdf.empty:
        print(
            f"    Error: empty shapefile for {basin}"
        )
        return

    # Cartopy PlateCarree expects longitude/latitude coordinates
    if boundary_gdf.crs is not None:
        boundary_gdf = boundary_gdf.to_crs(
            epsg=4326
        )

    (
        lon_min,
        lat_min,
        lon_max,
        lat_max,
    ) = boundary_gdf.total_bounds

    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min

    if lon_range <= 0 or lat_range <= 0:
        print(
            f"    Error: invalid spatial extent for "
            f"{basin}"
        )
        return

    # Add a small map margin
    lon_margin = lon_range * 0.02
    lat_margin = lat_range * 0.02

    plot_lon_min = lon_min - lon_margin
    plot_lon_max = lon_max + lon_margin
    plot_lat_min = lat_min - lat_margin
    plot_lat_max = lat_max + lat_margin

    plot_lon_range = plot_lon_max - plot_lon_min
    plot_lat_range = plot_lat_max - plot_lat_min

    aspect_ratio = (
        plot_lon_range / plot_lat_range
    )

    # --------------------------------------------------------
    # Set up figure
    # --------------------------------------------------------

    row_height = 5.0

    figure_width = max(
        row_height * aspect_ratio,
        6.0,
    )

    # Prevent extremely wide figures
    figure_width = min(
        figure_width,
        18.0,
    )

    figure_height = (
        row_height * len(TARGET_CROPS)
    )

    fig = plt.figure(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    grid_spec = gridspec.GridSpec(
        nrows=len(TARGET_CROPS),
        ncols=1,
        figure=fig,
        hspace=0.03,
        left=0.02,
        right=0.98,
        bottom=0.10,
        top=0.94,
    )

    fig.suptitle(
        f"{yield_mode} yield: {basin}",
        fontsize=30,
        fontweight="bold",
        y=0.975,
    )

    low_runoff = load_low_runoff_mask(
        basin
    )

    last_image = None
    available_crop_count = 0

    # --------------------------------------------------------
    # Plot the four crops
    # --------------------------------------------------------

    for crop_index, crop_name in enumerate(
        TARGET_CROPS
    ):

        ax = fig.add_subplot(
            grid_spec[crop_index, 0],
            projection=ccrs.PlateCarree(),
        )

        crop_file_names = CROP_MAP[
            crop_name
        ]

        yield_data = get_aggregated_yield(
            folder=mode_config["folder"],
            basin=basin,
            crop_list=crop_file_names,
            components=mode_config["components"],
        )

        if yield_data is not None:

            available_crop_count += 1

            last_image = yield_data.plot(
                ax=ax,
                transform=ccrs.PlateCarree(),
                cmap=yield_cmap,
                norm=yield_norm,
                add_colorbar=False,
                zorder=1,
            )

            low_runoff_for_plot = (
                prepare_low_runoff_mask(
                    low_runoff=low_runoff,
                    yield_data=yield_data,
                )
            )

            if low_runoff_for_plot is not None:

                low_runoff_for_plot.plot(
                    ax=ax,
                    transform=ccrs.PlateCarree(),
                    cmap=grey_cmap,
                    vmin=0,
                    vmax=1,
                    add_colorbar=False,
                    alpha=0.5,
                    zorder=2,
                )

        else:

            ax.text(
                0.5,
                0.5,
                "No data available",
                transform=ax.transAxes,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=22,
                color="black",
                zorder=5,
            )

            print(
                f"    No {yield_mode.lower()} "
                f"data found for {crop_name}"
            )

        # Basin outline
        boundary_gdf.boundary.plot(
            ax=ax,
            color="black",
            linewidth=0.8,
            zorder=3,
        )

        ax.set_extent(
            [
                plot_lon_min,
                plot_lon_max,
                plot_lat_min,
                plot_lat_max,
            ],
            crs=ccrs.PlateCarree(),
        )

        ax.axis("off")


        # Add the scale bar to the bottom crop panel only
        if crop_index == len(TARGET_CROPS) - 1:

            if plot_lon_range > 10:
                scale_length_km = 500
            elif plot_lon_range > 5:
                scale_length_km = 200
            else:
                scale_length_km = 100

            add_scale_bar(
                ax=ax,
                length_km=scale_length_km,
                basin=basin,
            )

    # --------------------------------------------------------
    # Add one shared colorbar
    # --------------------------------------------------------

    if last_image is not None:

        colorbar_axis = fig.add_axes(
            [
                0.18,
                0.045,
                0.64,
                0.018,
            ]
        )

        colorbar = fig.colorbar(
            last_image,
            cax=colorbar_axis,
            orientation="horizontal",
            extend="max",
        )

        colorbar.set_ticks([0, 2000, 4000, 6000, 8000, 10000])

        colorbar.ax.tick_params(
            labelsize=18,
            length=5,
        )

        colorbar.set_label(
            "Yield (kg ha$^{-1}$)",
            fontsize=24,
            labelpad=8,
        )

        colorbar.outline.set_visible(
            False
        )

    # --------------------------------------------------------
    # Save figure
    # --------------------------------------------------------

    output_filename = (
        f"{yield_mode}_Yield_{basin}.png"
    )

    output_path = os.path.join(
        FIG_OUT,
        output_filename,
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.05,
    )

    plt.close(fig)

    print(
        f"    Saved: {output_path}"
    )

    print(
        f"    Crop panels with data: "
        f"{available_crop_count}/"
        f"{len(TARGET_CROPS)}"
    )


# ============================================================
# 5. RUN SCRIPT
# ============================================================

def main():
    """
    Generate irrigated, rainfed, and potential-yield figures for all basins.
    """

    print("=" * 70)
    print("Generating basin yield figures")
    print("=" * 70)

    generated_figure_count = 0

    for basin in STUDY_AREAS:

        print(f"\nProcessing basin: {basin}")

        for yield_mode, mode_config in YIELD_MODES.items():

            plot_basin_yield(
                basin=basin,
                yield_mode=yield_mode,
                mode_config=mode_config,
            )

            generated_figure_count += 1

    print("\n" + "=" * 70)

    print(
        f"Finished. Requested figures: "
        f"{generated_figure_count}"
    )

    print(
        f"Output directory:\n{FIG_OUT}"
    )

    print("=" * 70)


if __name__ == "__main__":
    main()