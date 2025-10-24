# This code is used to check:
# 1. Where does the irrigation exceed the sustainable level
# 2. How much yield is contributed by these pixels

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

IrriDir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
OutputDir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs"
Basins = ["Indus", "Rhine", "LaPlata", "Yangtze"]
CropTypes = ["mainrice","secondrice", "maize", "winterwheat", "soybean"]
PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Part1"

start_year = 1986
end_year = 2015

# Pick colors for crops
crop_colors = {
    "mainrice": "#f19c9c",
    "secondrice": "#f19c9c",
    "maize": "#fcd562",
    "winterwheat": "#a2a1ef",
    "soybean": "#66c2a5"
}

for basin in Basins:
    print(f"=== Processing basin: {basin} ===")

    # Initialize lists to store averages for the bar chart
    avg_total_list = []
    avg_exceed_list = []
    crop_labels = []

    for crop in CropTypes:
        print(f"  > {crop}")

        # Load sustainable/unsustainable irrigation rate
        try:
            Unsus_Irri = xr.open_dataset(os.path.join(IrriDir, basin, "Unsus_Irrigation", f"{basin}_{crop}_monthly_Irri_Rate.nc")).sel(time=slice(str(start_year), str(end_year)))
            Sus_Irri = xr.open_dataset(os.path.join(IrriDir, basin, "Sus_Irrigation", f"{basin}_{crop}_monthly_Irri_Rate.nc")).sel(time=slice(str(start_year), str(end_year)))
        except FileNotFoundError:
            print(f"⚠️ Boundary files missing for {basin}, skipping basin.")
            continue
        Unsus_Irri_Rate = Unsus_Irri["Irrigation_Rate"]
        Sus_Irri_Rate = Sus_Irri["Irrigation_Rate"]

        # Sum to annual 
        Unsus_Irri_annual = Unsus_Irri_Rate.groupby("time.year").sum(dim="time", skipna=True)
        Sus_Irri_annual   = Sus_Irri_Rate.groupby("time.year").sum(dim="time", skipna=True)
        Unsus_Irri_annual = Unsus_Irri_annual.rename(year="time")
        Sus_Irri_annual   = Sus_Irri_annual.rename(year="time")

        # Load output
        output_nc_path_irrigated = os.path.join(OutputDir, f"Output_Unsus_Irrigation/{basin}_{crop}_annual.nc")
        output_nc_path_rainfed = os.path.join(OutputDir, f"Output_Rainfed/{basin}_{crop}_annual.nc")
        
        # Check path
        if not os.path.exists(output_nc_path_irrigated):
            print(f"    ⚠️ Missing output file: {output_nc_path_irrigated}")
            continue
        if not os.path.exists(output_nc_path_rainfed):
            print(f"    ⚠️ Missing output file: {output_nc_path_rainfed}")
            continue

        output_nc_irrigated = xr.open_dataset(output_nc_path_irrigated).sel(time=slice(start_year, end_year))
        output_nc_rainfed = xr.open_dataset(output_nc_path_rainfed).sel(time=slice(start_year, end_year))

        # Check variables
        required_vars = ["Storage"]
        if not all(var in output_nc_irrigated for var in required_vars):
            print(f"    ⚠️ Missing required variables in {output_nc_path_irrigated}, skipping crop.")
            continue
        if not all(var in output_nc_rainfed for var in required_vars):
            print(f"    ⚠️ Missing required variables in {output_nc_path_rainfed}, skipping crop.")
            continue

        # Load harvested area for irrigated 
        HA_path_irrigated = os.path.join(IrriDir, basin, "Mask", f"{crop}_HA_Irrigated.nc")
        if not os.path.exists(HA_path_irrigated):
            print(f"    ⚠️ Missing mask file: {HA_path_irrigated}, skipping crop.")
            continue
        HA_nc_irrigated = xr.open_dataset(HA_path_irrigated)
        HA_irrigated = HA_nc_irrigated["area_total"] if "area_total" in HA_nc_irrigated else list(HA_nc_irrigated.data_vars.values())[0]

        # Load harvested area for rainfed cropland
        HA_path_rainfed = os.path.join(IrriDir, basin, "Mask", f"{crop}_HA_rainfed.nc")
        if not os.path.exists(HA_path_rainfed):
            print(f"    ⚠️ Missing mask file: {HA_path_rainfed}, skipping crop.")
            continue
        HA_nc_rainfed = xr.open_dataset(HA_path_rainfed)
        HA_rainfed = HA_nc_rainfed["area_total"] if "area_total" in HA_nc_rainfed else list(HA_nc_rainfed.data_vars.values())[0]
        
        # Total harvest area: Irrigated + Rainfed
        HA = HA_rainfed.fillna(0) + HA_irrigated.fillna(0)

        # Check if any HA > 2500
        if (HA > 2500).sum() == 0:
            print(f" ⚠️ No data > 2500 in {crop} for {basin}, skipping crop.")
            continue

        valid_HA_irrigated = xr.where(HA > 2500, HA_irrigated, np.nan)
        valid_HA_rainfed = xr.where(HA > 2500, HA_rainfed, np.nan)

        # Pixels with unsustainable irrigation rate
        mask_unsus_irrigation = xr.where((Unsus_Irri_annual > Sus_Irri_annual), 1, 0)

        # Total production (kg → tons)
        total_production_tons = (((output_nc_rainfed["Storage"] * valid_HA_rainfed).sum(dim=["lat", "lon"], skipna=True)) +  ((output_nc_irrigated["Storage"] * valid_HA_irrigated).sum(dim=["lat", "lon"], skipna=True)))/ 1000
        exceed_production_tons = ((output_nc_irrigated["Storage"] * valid_HA_irrigated * mask_unsus_irrigation).sum(dim=["lat", "lon"], skipna=True)) / 1000

        # Store average annual values for bar chart
        avg_total_list.append(float(total_production_tons.mean()))
        avg_exceed_list.append(float(exceed_production_tons.mean()))
        crop_labels.append(crop)

        # Step 3: Plot yearly time series
        years = output_nc_rainfed["time"].values
        plt.figure(figsize=(8, 5))
        plt.plot(years, total_production_tons, label="Total Production", linewidth=2)
        plt.plot(years, exceed_production_tons, label="Unsustainble Irrigation Area Production", linewidth=2)
        plt.title(f"{basin} - {crop} (1986-2015)")
        plt.ylabel("Production (tons/year)")
        plt.xlabel("Year")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        fig_path = os.path.join(PlotDir, f"{basin}_{crop}_unsus_irri_production_1986-2015.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved yearly figure: {fig_path}")

    # Step 4: Create overlayed bar chart for the basin
    if len(avg_total_list) == 0:
        print(f"    ⚠️ No crops with data for {basin}, skipping bar chart.")
        continue

    x = np.arange(len(crop_labels))  # positions
    width = 0.6  # bar width

    plt.figure(figsize=(8, 5))
    # Light bars: total production
    plt.bar(x, avg_total_list, width, color=[crop_colors[c] for c in crop_labels], alpha=0.5, label="Total Production")
    # Dark bars: exceedance production on top (overlay)
    plt.bar(x, avg_exceed_list, width, color=[crop_colors[c] for c in crop_labels], alpha=1.0, label="Unsustainble irrigated area")

    plt.xticks(x, crop_labels)
    plt.ylabel("Average Annual Production (tons/year)")
    plt.title(f"{basin} - Average Annual Production (1986-2015)")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    fig_path_bar = os.path.join(PlotDir, f"{basin}_Irri_Sum_barchart.png")
    plt.savefig(fig_path_bar, dpi=300)
    plt.close()
    print(f"    ✅ Saved overlayed summary bar chart: {fig_path_bar}")

print("=== All processing complete ===")