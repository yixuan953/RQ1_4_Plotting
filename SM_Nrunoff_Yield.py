import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


# ============================================================
# CONFIGURATION
# ============================================================

SCENARIOS = {
    "Current": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline",
    "Water": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_2_Sus_Irrigation",
    "Water + N": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_N_50",
    "Water + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_P_50",
    "Water + N + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_NP_50",
}

SUMMARY_DIR = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
)

BOUNDARY_DIR = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "4_RQ1_Analysis_Results/V5_Statistics/1_Boundary"
)

OUT_DIR = (
    "/lustre/nobackup/WUR/ESG/zhou111/"
    "4_RQ1_Analysis_Results/V5_Plots/SM_Fig/"
    "Annual_Variation/Nrunoff_Yield"
)

os.makedirs(OUT_DIR, exist_ok=True)

STUDY_AREAS = ["Rhine", "LaPlata", "Indus", "Yangtze"]
CROPS = ["winterwheat", "mainrice", "maize", "soybean"]

# Basin-specific reporting groups.
# For Yangtze, "Rice" combines mainrice and secondrice.
CROP_GROUPS = {
    "Yangtze": {
        "winterwheat": ["winterwheat"],
        "Rice": ["mainrice", "secondrice"],
        "maize": ["maize"],
        "soybean": ["soybean"],
    }
}

COLORS = ["#7A7878", "#458FBA", "#C06497", "#C89119", "#06956F"]

START_YEAR = 2009
END_YEAR = 2020
MIN_TOTAL_HA = 2500

# Preserve the +1 year shift from your original script.
# Set to False if the CSV year should be plotted directly.
SHIFT_YEAR_BY_ONE = True

# Same boundary scaling used in your original P script.
N_BOUNDARY_FACTOR = 0.70


# ============================================================
# POSSIBLE COLUMN NAMES
# Add your exact CSV names here when necessary.
# ============================================================

N_SURFACE_CANDIDATES = [
    "N_surf",
    "N_surface",
    "N_runoff_surf",
    "N_runoff_surface",
    "SurfaceN",
]

N_SUBSURFACE_CANDIDATES = [
    "N_sub",
    "N_subsurface",
    "N_runoff_sub",
    "N_runoff_subsurface",
    "SubsurfaceN",
]

IRRIGATED_YIELD_CANDIDATES = [
    "Storage",
    "Yield_Irrigated",
    "Irrigated_Yield",
    "Yield_irri",
    "Yield_Irri",
    "yield_irrigated",
    "yield_irri",
    "Avg_Yield_Irrigated",
]

RAINFED_YIELD_CANDIDATES = [
    "Storage",
    "Yield_Rainfed",
    "Rainfed_Yield",
    "Yield_rf",
    "Yield_RF",
    "yield_rainfed",
    "yield_rf",
    "Avg_Yield_Rainfed",
]


# ============================================================
# HELPERS
# ============================================================

def find_column(df, candidates, description):
    """Find the first matching column, ignoring case."""
    for name in candidates:
        if name in df.columns:
            return name

    lookup = {
        str(column).strip().lower(): column
        for column in df.columns
    }

    for name in candidates:
        match = lookup.get(name.strip().lower())
        if match is not None:
            return match

    warnings.warn(
        f"Could not find {description}. Tried {candidates}. "
        f"Available columns: {list(df.columns)}"
    )
    return None


def to_numeric(df, columns):
    for column in columns:
        if column is not None and column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")


def get_crop_groups(basin):
    """
    Return reporting-crop names and the source crop files used for each one.

    Example for Yangtze:
        "Rice" -> ["mainrice", "secondrice"]
    """
    if basin in CROP_GROUPS:
        return CROP_GROUPS[basin]

    return {
        crop: [crop]
        for crop in CROPS
    }


def combine_result_tables(result_tables):
    """
    Sum annual values from multiple component crops.

    Used to combine Yangtze mainrice and secondrice into one Rice series.
    """
    valid_tables = [
        table
        for table in result_tables
        if table is not None and not table.empty
    ]

    if not valid_tables:
        return pd.DataFrame(columns=["Year", "Value"])

    combined = pd.concat(
        valid_tables,
        ignore_index=True,
    )

    combined = (
        combined
        .groupby("Year", as_index=False)["Value"]
        .sum()
        .sort_values("Year")
    )

    return combined


def load_harvested_area(basin, crop):
    """Load valid irrigated and rainfed harvested areas."""
    summary_file = os.path.join(
        SUMMARY_DIR,
        f"{basin}_{crop}_summary.nc",
    )

    if not os.path.exists(summary_file):
        print(f"  Missing summary file: {summary_file}")
        return None

    with xr.open_dataset(summary_file) as ds:
        required = ["Total_HA", "Irrigated_HA", "Rainfed_HA"]
        missing = [name for name in required if name not in ds]

        if missing:
            print(f"  Missing variables in summary file: {missing}")
            return None

        valid = ds["Total_HA"] > MIN_TOTAL_HA

        area_ds = xr.Dataset({
            "HA_irri": ds["Irrigated_HA"].where(valid),
            "HA_rf": ds["Rainfed_HA"].where(valid),
        })

        area_df = area_ds.to_dataframe().reset_index()

    if not {"lat", "lon"}.issubset(area_df.columns):
        print(
            f"  Expected lat/lon coordinates, found: "
            f"{list(area_df.columns)}"
        )
        return None

    area_df = area_df.dropna(
        subset=["HA_irri", "HA_rf"],
        how="all",
    )

    area_df["HA_irri"] = area_df["HA_irri"].fillna(0)
    area_df["HA_rf"] = area_df["HA_rf"].fillna(0)

    return area_df


def load_n_boundary(basin, crops):
    """
    Read and sum crop-specific N boundaries.

    Parameters
    ----------
    basin : str
        Basin name.
    crops : str or list[str]
        One crop name or multiple component crops.

    For Yangtze Rice, pass:
        ["mainrice", "secondrice"]

    The returned safe boundary is:
        sum(component N boundaries) * N_BOUNDARY_FACTOR
    """
    if isinstance(crops, str):
        crops = [crops]

    boundary_file = os.path.join(
        BOUNDARY_DIR,
        f"{basin}_crop-specific_boundaries.csv",
    )

    if not os.path.exists(boundary_file):
        print(f"  Missing boundary file: {boundary_file}")
        return None

    bdf = pd.read_csv(boundary_file)

    if bdf.empty:
        return None

    crop_column = bdf.columns[0]

    boundary_column = find_column(
        bdf,
        [
            "N [ktons]",
            "N [kton]",
            "N[ktons]",
            "N_boundary [ktons]",
        ],
        "N boundary column",
    )

    if boundary_column is None:
        return None

    normalized_names = (
        bdf[crop_column]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    boundary_values = []

    for crop in crops:
        row = bdf[
            normalized_names == crop.strip().lower()
        ]

        if row.empty:
            print(
                f"  No N-boundary row for '{crop}' "
                f"in {boundary_file}"
            )
            continue

        value = pd.to_numeric(
            row.iloc[0][boundary_column],
            errors="coerce",
        )

        if pd.notna(value):
            boundary_values.append(float(value))

    if not boundary_values:
        return None

    return sum(boundary_values) * N_BOUNDARY_FACTOR


def load_and_merge(csv_file, area_df):
    """Load scenario CSV and merge it to the harvested-area grid."""
    df = pd.read_csv(csv_file)

    required = ["Year", "Lat", "Lon"]
    missing = [name for name in required if name not in df]

    if missing:
        print(f"  Missing CSV columns {missing}: {csv_file}")
        return None

    to_numeric(df, required)

    df = df[
        df["Year"].between(START_YEAR, END_YEAR)
    ].copy()

    merged = df.merge(
        area_df,
        left_on=["Lat", "Lon"],
        right_on=["lat", "lon"],
        how="inner",
    )

    if merged.empty:
        print(f"  No matching Lat/Lon cells: {csv_file}")
        return None

    return merged


def calculate_n_runoff(merged):
    """
    Total annual N runoff:
      sum((N_surf + N_sub) * irrigated HA)
      + sum((N_surf + N_sub) * rainfed HA)

    The 1e-6 conversion is retained from the original P calculation.
    """
    surf_col = find_column(
        merged,
        N_SURFACE_CANDIDATES,
        "surface N runoff",
    )
    sub_col = find_column(
        merged,
        N_SUBSURFACE_CANDIDATES,
        "subsurface N runoff",
    )

    if surf_col is None or sub_col is None:
        return pd.DataFrame(columns=["Year", "Value"])

    to_numeric(
        merged,
        [surf_col, sub_col, "HA_irri", "HA_rf"],
    )

    results = []

    for year, group in merged.groupby("Year"):
        n_per_ha = (
            group[surf_col].fillna(0)
            + group[sub_col].fillna(0)
        )

        total = (
            (n_per_ha * group["HA_irri"].fillna(0)).sum()
            + (n_per_ha * group["HA_rf"].fillna(0)).sum()
        ) * 1e-6

        results.append((year, total))

    return pd.DataFrame(
        results,
        columns=["Year", "Value"],
    ).sort_values("Year")


def calculate_combined_yield(merged):
    """
    Calculate total basin-scale crop production from irrigated and
    rainfed systems using the CSV variable "Storage".

    Basin crop production =
        sum(irrigated yield * irrigated harvested area)
        +
        sum(rainfed yield * rainfed harvested area)

    Assuming:
        yield is in kg/ha
        harvested area is in ha

    The product is kg/year and is converted to Mt/year using:

        1 Mt = 1e9 kg
    """
    yield_irri_col = find_column(
        merged,
        IRRIGATED_YIELD_CANDIDATES,
        "irrigated yield",
    )
    yield_rf_col = find_column(
        merged,
        RAINFED_YIELD_CANDIDATES,
        "rainfed yield",
    )

    if yield_irri_col is None or yield_rf_col is None:
        return pd.DataFrame(columns=["Year", "Value"])

    to_numeric(
        merged,
        [
            yield_irri_col,
            yield_rf_col,
            "HA_irri",
            "HA_rf",
        ],
    )

    results = []

    for year, group in merged.groupby("Year"):
        irrigated_production_kg = (
            group[yield_irri_col].fillna(0)
            * group["HA_irri"].fillna(0)
        ).sum()

        rainfed_production_kg = (
            group[yield_rf_col].fillna(0)
            * group["HA_rf"].fillna(0)
        ).sum()

        total_production_mt = (
            irrigated_production_kg
            + rainfed_production_kg
        ) * 1e-9

        results.append((year, total_production_mt))

    return pd.DataFrame(
        results,
        columns=["Year", "Value"],
    ).sort_values("Year")


def x_values(values):
    if SHIFT_YEAR_BY_ONE:
        return values["Year"] + 1
    return values["Year"]


def set_nice_yaxis(ax, n_ticks=5):
    """
    Set exactly five evenly spaced y-axis tick labels using a rounded
    interval based on the plotted data range.
    """
    ymin, ymax = ax.get_ylim()

    if not np.isfinite(ymin) or not np.isfinite(ymax):
        return

    if np.isclose(ymin, ymax):
        padding = abs(ymin) * 0.1 if ymin != 0 else 1.0
        ymin -= padding
        ymax += padding

    raw_step = (ymax - ymin) / (n_ticks - 1)

    if raw_step <= 0:
        return

    magnitude = 10 ** np.floor(np.log10(raw_step))
    normalized_step = raw_step / magnitude

    if normalized_step <= 1:
        nice_normalized_step = 1
    elif normalized_step <= 2:
        nice_normalized_step = 2
    elif normalized_step <= 2.5:
        nice_normalized_step = 2.5
    elif normalized_step <= 5:
        nice_normalized_step = 5
    else:
        nice_normalized_step = 10

    nice_step = nice_normalized_step * magnitude

    nice_min = np.floor(ymin / nice_step) * nice_step
    nice_max = nice_min + nice_step * (n_ticks - 1)

    # Shift the five-tick range upward until it includes the data maximum.
    while nice_max < ymax:
        nice_min += nice_step
        nice_max += nice_step

    ticks = nice_min + nice_step * np.arange(n_ticks)

    ax.set_ylim(ticks[0], ticks[-1])
    ax.set_yticks(ticks)


def style_axis(ax, ylabel):
    ax.set_xlabel("Year", fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.grid(True, linestyle="--", alpha=0.5)


def plot_n_runoff(basin, crop, results, boundary):
    """Save a separate N-runoff figure."""
    fig, ax = plt.subplots(figsize=(7, 4.5))

    if boundary is not None:
        ax.axhline(
            boundary,
            color="red",
            linestyle="--",
            linewidth=3,
            label="Safe boundary",
        )

    for index, scenario in enumerate(SCENARIOS):
        values = results.get(scenario)

        if values is None or values.empty:
            continue

        ax.plot(
            x_values(values),
            values["Value"],
            color=COLORS[index],
            linewidth=2.5,
            label=scenario,
        )

    style_axis(ax, "N runoff [kton N/yr]")
    set_nice_yaxis(ax)
    fig.tight_layout()

    output_file = os.path.join(
        OUT_DIR,
        f"{basin}_{crop}_N_runoff.png",
    )

    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_file}")


def plot_combined_yield(basin, crop, results):
    """Save basin-scale irrigated + rainfed crop production."""
    fig, ax = plt.subplots(figsize=(7, 4.5))

    for index, scenario in enumerate(SCENARIOS):
        values = results.get(scenario)

        if values is None or values.empty:
            continue

        ax.plot(
            x_values(values),
            values["Value"],
            color=COLORS[index],
            linewidth=2.5,
            label=scenario,
        )

    style_axis(ax, "Crop yield [Mt yr$^{-1}$]")
    set_nice_yaxis(ax)
    fig.tight_layout()

    output_file = os.path.join(
        OUT_DIR,
        f"{basin}_{crop}_combined_yield.png",
    )

    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_file}")


# ============================================================
# MAIN LOOP
# ============================================================

def main():
    """
    Generate:
      1. N-runoff figure with boundary.
      2. Basin-scale irrigated + rainfed crop-production figure.

    For Yangtze, the reporting crop "Rice" is the annual sum of:
      - mainrice
      - secondrice

    Its safe boundary is the sum of the mainrice and secondrice rows.
    """
    for basin in STUDY_AREAS:

        crop_groups = get_crop_groups(basin)

        for report_crop, source_crops in crop_groups.items():

            print(
                f"\nProcessing {basin} - {report_crop} "
                f"from {source_crops}"
            )

            # Sum the component-crop boundaries.
            # Yangtze Rice therefore uses mainrice + secondrice.
            boundary = load_n_boundary(
                basin,
                source_crops,
            )

            runoff_results = {}
            yield_results = {}

            for scenario, scenario_dir in SCENARIOS.items():

                scenario_runoff_tables = []
                scenario_yield_tables = []

                for source_crop in source_crops:

                    area_df = load_harvested_area(
                        basin,
                        source_crop,
                    )

                    if area_df is None or area_df.empty:
                        print(
                            f"  No harvested-area data for "
                            f"{source_crop}"
                        )
                        continue

                    csv_file = os.path.join(
                        scenario_dir,
                        f"{basin}_{source_crop}_annual.csv",
                    )

                    if not os.path.exists(csv_file):
                        print(
                            f"  Missing scenario CSV: "
                            f"{csv_file}"
                        )
                        continue

                    print(
                        f"  Scenario: {scenario}; "
                        f"source crop: {source_crop}"
                    )

                    merged = load_and_merge(
                        csv_file,
                        area_df,
                    )

                    if merged is None or merged.empty:
                        continue

                    scenario_runoff_tables.append(
                        calculate_n_runoff(
                            merged.copy()
                        )
                    )

                    scenario_yield_tables.append(
                        calculate_combined_yield(
                            merged.copy()
                        )
                    )

                runoff_results[scenario] = (
                    combine_result_tables(
                        scenario_runoff_tables
                    )
                )

                yield_results[scenario] = (
                    combine_result_tables(
                        scenario_yield_tables
                    )
                )

            plot_n_runoff(
                basin,
                report_crop,
                runoff_results,
                boundary,
            )

            plot_combined_yield(
                basin,
                report_crop,
                yield_results,
            )


if __name__ == "__main__":
    main()