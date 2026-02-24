import pandas as pd
import xarray as xr
import numpy as np
import os

# study areas and crop types
studyareas = ["Indus", "Rhine", "LaPlata", "Yangtze"]
croptypes = ["mainrice", "secondrice", "maize", "soybean", "winterwheat"]


# Irrigated paths
csv_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Inc_prop"
range_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
out_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Inc_prop"

# Rainfed paths
# csv_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed/Inc_prop"
# range_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
# out_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed/Inc_prop"


for studyarea in studyareas:
    # load reference grid
    range_file = f"{range_dir}/{studyarea}/range.nc"
    if not os.path.exists(range_file):
        print(f"!!! range.nc missing for {studyarea}, skipping...")
        continue

    ref = xr.open_dataset(range_file)
    lat = ref["lat"].values
    lon = ref["lon"].values

    # precompute index maps
    lat_index = {v: i for i, v in enumerate(lat)}
    lon_index = {v: j for j, v in enumerate(lon)}

    for croptype in croptypes:
        csv_file = f"{csv_dir}/{studyarea}_{croptype}_annual.csv"
        output_file = f"{out_dir}/{studyarea}_{croptype}_annual.nc"

        if not os.path.exists(csv_file):
            print(f" No CSV for {studyarea} - {croptype}, skipping...")
            continue

        print(f"Processing {studyarea} - {croptype} ...")

        # read CSV
        df = pd.read_csv(csv_file)

        years = sorted(df["Year"].unique())
        year_index = {v: k for k, v in enumerate(years)}

        # init arrays (time, lat, lon)
        shape = (len(years), len(lat), len(lon))
        storage_arr = np.full(shape, np.nan)
        growthday_arr = np.full(shape, np.nan)
        harvestday_arr = np.full(shape, np.nan)
        nup_arr = np.full(shape, np.nan)
        pup_arr = np.full(shape, np.nan)
        nloss_arr = np.full(shape, np.nan)
        ploss_arr = np.full(shape, np.nan)
        ngrain_arr = np.full(shape, np.nan)
        pgrain_arr = np.full(shape, np.nan)

        # fill arrays
        for _, row in df.iterrows():
            i = lat_index.get(row["Lat"])
            j = lon_index.get(row["Lon"])
            t = year_index.get(row["Year"])
            if i is not None and j is not None and t is not None:
                storage_arr[t, i, j] = row["Storage"]
                growthday_arr[t, i, j] = row["GrowthDay"]
                harvestday_arr[t, i, j] = row["Day"]
                nup_arr[t, i, j] = row["N_uptake"]
                pup_arr[t, i, j] = row["P_uptake"]
                nloss_arr[t, i, j] = row["N_surf"] + row["N_sub"] 
                ploss_arr[t, i, j] = row["P_surf"] + row["P_sub"] 
                ngrain_arr[t, i, j] = row["N_grain"]
                pgrain_arr[t, i, j] = row["P_grain"]

        # make Dataset
        ds = xr.Dataset(
            {
                "Yield": (("year", "lat", "lon"), storage_arr),
                "GrowthDay": (("year", "lat", "lon"), growthday_arr),
                "HarvestDay": (("year", "lat", "lon"), harvestday_arr),
                "N_Uptake": (("year", "lat", "lon"), nup_arr),
                "P_Uptake": (("year", "lat", "lon"), pup_arr),
                "N_Runoff": (("year", "lat", "lon"), nloss_arr),
                "P_Runoff": (("year", "lat", "lon"), ploss_arr),
                "N_grain": (("year", "lat", "lon"), ngrain_arr),
                "P_grain": (("year", "lat", "lon"), pgrain_arr),

            },
            coords={"year": years, "lat": lat, "lon": lon}
        )

        # attributes
        ds.attrs["description"] = f"Annual yield for {croptype} in {studyarea}"
        for var in ds.data_vars:
            ds[var].attrs["_FillValue"] = 0.0

        # replace NaN with FillValue
        ds = ds.fillna(0.0)

        # save
        ds.to_netcdf(output_file)
        print(f" Saved {output_file}")
