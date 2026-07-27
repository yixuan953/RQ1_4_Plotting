# This script is used to calculate how many people are living in the boundaries exceeded pixels


import xarray as xr
import rasterio
import numpy as np
import pandas as pd
import os

# --------------------------------------------------
# Files
# --------------------------------------------------
pop_file = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_min.tif"

studyareas = ["LaPlata", "Yangtze", "Indus", "Rhine"]

out_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/Extended_Stat/3_Boundaries_Exceeded_pop/1_Baseline"
# out_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/Extended_Stat/3_Boundaries_Exceeded_pop/2_RespectNPWater"
os.makedirs(out_dir, exist_ok=True)

# All classes
classes = [111, 112, 121, 122, 211, 212, 221, 222]

# Group: everything except 111
exceed_classes = [112, 121, 122, 211, 212, 221, 222]

# --------------------------------------------------
# Open population raster once
# --------------------------------------------------
with rasterio.open(pop_file) as src:
    pop_data = src.read(1)
    pop_transform = src.transform
    pop_nodata = src.nodata

    if pop_nodata is not None:
        pop_data = np.where(pop_data == pop_nodata, np.nan, pop_data)

# --------------------------------------------------
# Helper: sample population to basin grid
# --------------------------------------------------
def sample_population(lon, lat):

    cols, rows = np.meshgrid(lon, lat)

    px, py = ~pop_transform * (cols, rows)

    px = px.astype(int)
    py = py.astype(int)

    px = np.clip(px, 0, pop_data.shape[1] - 1)
    py = np.clip(py, 0, pop_data.shape[0] - 1)

    return pop_data[py, px]

# --------------------------------------------------
# Loop basins
# --------------------------------------------------
for basin in studyareas:

    print(f"\nProcessing {basin} ...")

    nc_file = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline/{basin}_boundary_check.nc"

    out_csv = f"{out_dir}/{basin}_population_by_boundary_class.csv"

    # --------------------------------------------------
    # Read NetCDF
    # --------------------------------------------------
    ds = xr.open_dataset(nc_file)

    varname = [v for v in ds.data_vars if v != "spatial_ref"][0]

    boundary = ds[varname].values
    lon = ds["x"].values
    lat = ds["y"].values

    print("Boundary shape:", boundary.shape)

    # --------------------------------------------------
    # Sample population
    # --------------------------------------------------
    pop_local = sample_population(lon, lat)

    print("Sampled population shape:", pop_local.shape)

    # --------------------------------------------------
    # Class-wise population
    # --------------------------------------------------
    results = []

    for c in classes:

        mask = boundary == c

        pop_sum = np.nansum(pop_local[mask])
        n_cells = np.sum(mask)

        results.append({
            "Class": c,
            "Cells": int(n_cells),
            "Population": float(pop_sum)
        })

    df = pd.DataFrame(results)

    # --------------------------------------------------
    # Grouped summaries
    # --------------------------------------------------
    total_population = np.nansum(pop_local[~np.isnan(boundary)])
    exceed_population = np.nansum(pop_local[np.isin(boundary, exceed_classes)])

    summary = pd.DataFrame([{
        "Class": "ALL_VALID",
        "Cells": int(np.sum(~np.isnan(boundary))),
        "Population": float(total_population)
    },
    {
        "Class": "EXCEED_CLASSES_112_222",
        "Cells": int(np.sum(np.isin(boundary, exceed_classes))),
        "Population": float(exceed_population)
    }])

    final_df = pd.concat([df, summary], ignore_index=True)

    print(final_df)

    final_df.to_csv(out_csv, index=False)

    print(f"Saved: {out_csv}")