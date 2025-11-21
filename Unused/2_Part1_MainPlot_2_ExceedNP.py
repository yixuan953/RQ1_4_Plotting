import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# === Paths ===
BaseDir = "/lustre/nobackup/WUR/ESG/zhou111"
ExceedDir = f"{BaseDir}/3_RQ1_Model_Outputs/Test_CriticalNP/Method1"
YieldDir_Irrig = f"{BaseDir}/3_RQ1_Model_Outputs/Output_Unsus_Irrigation"
YieldDir_Rain = f"{BaseDir}/3_RQ1_Model_Outputs/Output_Rainfed"
DataDir = f"{BaseDir}/2_RQ1_Data/2_StudyArea"
PlotDir = f"{BaseDir}/4_RQ1_Analysis_Results/Part1"

Basins = ["Indus", "Rhine", "LaPlata", "Yangtze"]
CropTypes = ["mainrice", "maize", "winterwheat", "soybean", "secondrice"]

start_year = 2005
end_year = 2015

# === Crop colors ===
crop_colors = {
    "mainrice": "#f19c9c",
    "secondrice": "#f19c9c",
    "maize": "#fcd562",
    "winterwheat": "#a2a1ef",
    "soybean": "#66c2a5"
}

for basin in Basins:
    print(f"\n=== Processing {basin} ===")

    # --- Load exceedance (polluted) pixel mask ---
    exceed_path = os.path.join(ExceedDir, f"{basin}_exceed_pixels.nc")
    if not os.path.exists(exceed_path):
        print(f"⚠️ Missing exceed file for {basin}, skipping.")
        continue
    exceed_mask = xr.open_dataset(exceed_path)
    exceed_mask = list(exceed_mask.data_vars.values())[0]  # assume 1 var
    exceed_mask = exceed_mask.where(exceed_mask == 1, 0)   # ensure binary mask

    # --- Initialize plot data ---
    avg_total_list = []
    avg_polluted_list = []
    crop_labels = []

    for crop in CropTypes:
        print(f"  > {crop}")

        # === Load yield (irrigated and rainfed) ===
        yield_irrig_path = os.path.join(YieldDir_Irrig, f"{basin}_{crop}_annual.nc")
        yield_rain_path = os.path.join(YieldDir_Rain, f"{basin}_{crop}_annual.nc")
        if not (os.path.exists(yield_irrig_path) and os.path.exists(yield_rain_path)):
            print(f"    ⚠️ Missing yield files for {crop}, skipping.")
            continue

        yield_irrig = xr.open_dataset(yield_irrig_path).sel(time=slice(start_year, end_year))
        yield_rain = xr.open_dataset(yield_rain_path).sel(time=slice(start_year, end_year))

        if "Storage" not in yield_irrig or "Storage" not in yield_rain:
            print(f"    ⚠️ Missing 'Storage' variable in yield files, skipping.")
            continue

        Y_irrig = yield_irrig["Storage"]
        Y_rain = yield_rain["Storage"]

        # === Load harvested areas ===
        HA_rainfed = xr.open_dataset(os.path.join(DataDir, basin, "Mask", f"{crop}_HA_rainfed.nc"))
        HA_rainfed = HA_rainfed["area_total"] if "area_total" in HA_rainfed else list(HA_rainfed.data_vars.values())[0]

        HA_irrigated = xr.open_dataset(os.path.join(DataDir, basin, "Mask", f"{crop}_HA_Irrigated.nc"))
        HA_irrigated = HA_irrigated["area_total"] if "area_total" in HA_irrigated else list(HA_irrigated.data_vars.values())[0]

        # === Compute production (tons) ===
        prod_rain = (Y_rain * HA_rainfed).sum(dim=["lat", "lon"], skipna=True) / 1000
        prod_irrig = (Y_irrig * HA_irrigated).sum(dim=["lat", "lon"], skipna=True) / 1000
        total_prod = prod_rain + prod_irrig

        # === Polluted production (tons) ===
        prod_rain_polluted = (Y_rain * HA_rainfed * exceed_mask).sum(dim=["lat", "lon"], skipna=True) / 1000
        prod_irrig_polluted = (Y_irrig * HA_irrigated * exceed_mask).sum(dim=["lat", "lon"], skipna=True) / 1000
        polluted_prod = prod_rain_polluted + prod_irrig_polluted

        # === Average for 2005–2015 ===
        avg_total = float(total_prod.mean())
        avg_polluted = float(polluted_prod.mean())

        avg_total_list.append(avg_total)
        avg_polluted_list.append(avg_polluted)
        crop_labels.append(crop)

    # --- Plot bar chart ---
    if not avg_total_list:
        print(f"    ⚠️ No data to plot for {basin}.")
        continue

    x = np.arange(len(crop_labels))
    width = 0.6

    plt.figure(figsize=(8, 5))
    plt.bar(x, avg_total_list, width, color=[crop_colors[c] for c in crop_labels], alpha=0.5, label="Total Production")
    plt.bar(x, avg_polluted_list, width, color=[crop_colors[c] for c in crop_labels], alpha=1.0, label="Polluted Pixels Production")

    plt.xticks(x, crop_labels)
    plt.ylabel("Average Annual Production (tons/year)")
    plt.title(f"{basin} — 2005–2015 Average Production")
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()

    out_path = os.path.join(PlotDir, f"{basin}_production_polluted_vs_total_2005-2015.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"✅ Saved bar chart: {out_path}")

print("\n=== All basins processed ===")