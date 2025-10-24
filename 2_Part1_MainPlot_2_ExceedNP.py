# This code is used to:
# 1. Get the pixels where the N, P boudnaries have been exceeded for each year in 1986 - 2015
# 2. Calculate how much crop production was contributed by these pixels (by crop)
# 3. Compare the above-mentioned productions with total productions

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

BoundaryDir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
OutputDir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs"
Basins = ["Indus", "Rhine", "LaPlata", "Yangtze"]
CropTypes = ["mainrice", "maize", "winterwheat", "soybean", "secondrice"]
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

    # Load boundary datasets and subset to 1986-2015
    try:
        N_leach_boundary_nc = xr.open_dataset(os.path.join(BoundaryDir, basin, "Hydro", f"{basin}_N_critical_leaching_annual.nc")).sel(time=slice(str(start_year), str(end_year)))
        P_leach_boundary_nc = xr.open_dataset(os.path.join(BoundaryDir, basin, "Hydro", f"{basin}_P_critical_leaching_annual.nc")).sel(time=slice(str(start_year), str(end_year)))
        N_runoff_boundary_nc = xr.open_dataset(os.path.join(BoundaryDir, basin, "Hydro", f"{basin}_N_critical_runoff_annual.nc")).sel(time=slice(str(start_year), str(end_year)))
        P_runoff_boundary_nc = xr.open_dataset(os.path.join(BoundaryDir, basin, "Hydro", f"{basin}_P_critical_runoff_annual.nc")).sel(time=slice(str(start_year), str(end_year)))
    except FileNotFoundError:
        print(f"⚠️ Boundary files missing for {basin}, skipping basin.")
        continue

    # Apply concentration multipliers
    N_leach_boundary = N_leach_boundary_nc["OUT_BASEFLOW"]
    P_leach_boundary = P_leach_boundary_nc["OUT_BASEFLOW"]
    N_runoff_boundary = N_runoff_boundary_nc["OUT_RUNOFF"]
    P_runoff_boundary = P_runoff_boundary_nc["OUT_RUNOFF"]

    # Convert boundary times to integer years to align with output_nc
    N_leach_boundary = N_leach_boundary.assign_coords(time=N_leach_boundary["time.year"])
    P_leach_boundary = P_leach_boundary.assign_coords(time=P_leach_boundary["time.year"])
    N_runoff_boundary = N_runoff_boundary.assign_coords(time=N_runoff_boundary["time.year"])
    P_runoff_boundary = P_runoff_boundary.assign_coords(time=P_runoff_boundary["time.year"])

    for crop in CropTypes:
        print(f"  > {crop}")
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
        required_vars = ["N_leach", "N_surf", "N_sub", "P_leach", "P_surf", "P_sub", "Storage"]
        if not all(var in output_nc_irrigated for var in required_vars):
            print(f"    ⚠️ Missing required variables in {output_nc_path_irrigated}, skipping crop.")
            continue
        if not all(var in output_nc_rainfed for var in required_vars):
            print(f"    ⚠️ Missing required variables in {output_nc_path_rainfed}, skipping crop.")
            continue

        # Compute N, P losses for irrigated cropland
        N_leach_irrigated = output_nc_irrigated["N_leach"]
        N_runoff_irrigated = output_nc_irrigated["N_surf"].fillna(0) + output_nc_irrigated["N_sub"].fillna(0)
        P_leach_irrigated = output_nc_irrigated["P_leach"]
        P_runoff_irrigated = output_nc_irrigated["P_surf"].fillna(0) + output_nc_irrigated["P_sub"].fillna(0)
        # Load harvested area for irrigated 
        HA_path_irrigated = os.path.join(BoundaryDir, basin, "Mask", f"{crop}_HA_Irrigated.nc")
        if not os.path.exists(HA_path_irrigated):
            print(f"    ⚠️ Missing mask file: {HA_path_irrigated}, skipping crop.")
            continue
        HA_nc_irrigated = xr.open_dataset(HA_path_irrigated)
        HA_irrigated = HA_nc_irrigated["area_total"] if "area_total" in HA_nc_irrigated else list(HA_nc_irrigated.data_vars.values())[0]

        # Compute N, P losses for rainfed cropland
        N_leach_rainfed = output_nc_rainfed["N_leach"]
        N_runoff_rainfed = output_nc_rainfed["N_surf"].fillna(0) + output_nc_rainfed["N_sub"].fillna(0)
        P_leach_rainfed = output_nc_rainfed["P_leach"]
        P_runoff_rainfed = output_nc_rainfed["P_surf"].fillna(0) + output_nc_rainfed["P_sub"].fillna(0)
        # Load harvested area for irrigated 
        HA_path_rainfed = os.path.join(BoundaryDir, basin, "Mask", f"{crop}_HA_rainfed.nc")
        if not os.path.exists(HA_path_rainfed):
            print(f"    ⚠️ Missing mask file: {HA_path_rainfed}, skipping crop.")
            continue
        HA_nc_rainfed = xr.open_dataset(HA_path_rainfed)
        HA_rainfed = HA_nc_rainfed["area_total"] if "area_total" in HA_nc_rainfed else list(HA_nc_rainfed.data_vars.values())[0]

        HA = HA_rainfed.fillna(0) + HA_irrigated.fillna(0)

        # Check if any HA > 2500
        if (HA > 2500).sum() == 0:
            print(f" ⚠️ No data > 2500 in {crop} for {basin}, skipping crop.")
            continue

        valid_HA_irrigated = xr.where(HA > 2500, HA_irrigated, np.nan)
        valid_HA_rainfed = xr.where(HA > 2500, HA_rainfed, np.nan)

        # Exceedance mask (yearly)
        mask_data_irrigated = xr.where(
            (N_leach_irrigated > N_leach_boundary), # |
            # (P_leach_irrigated > P_leach_boundary) |
            # (N_runoff_irrigated > N_runoff_boundary) |
            # (P_runoff_irrigated > P_runoff_boundary),
            1, 0
        )

        mask_data_rainfed = xr.where(
            (N_leach_rainfed > N_leach_boundary), # |
            # (P_leach_rainfed > P_leach_boundary) |
            # (N_runoff_rainfed > N_runoff_boundary) |
            # (P_runoff_rainfed > P_runoff_boundary),
            1, 0
        )

        # Total production (kg → tons)
        total_production_tons = (((output_nc_rainfed["Storage"] * valid_HA_rainfed).sum(dim=["lat", "lon"], skipna=True)) +  ((output_nc_irrigated["Storage"] * valid_HA_irrigated).sum(dim=["lat", "lon"], skipna=True)))/ 1000
        exceed_production_tons = (((output_nc_rainfed["Storage"] * valid_HA_rainfed * mask_data_rainfed).sum(dim=["lat", "lon"], skipna=True)) + ((output_nc_irrigated["Storage"] * valid_HA_irrigated * mask_data_rainfed).sum(dim=["lat", "lon"], skipna=True))) / 1000

        # Store average annual values for bar chart
        avg_total_list.append(float(total_production_tons.mean()))
        avg_exceed_list.append(float(exceed_production_tons.mean()))
        crop_labels.append(crop)

        # Step 3: Plot yearly time series
        years = output_nc_rainfed["time"].values
        plt.figure(figsize=(8, 5))
        plt.plot(years, total_production_tons, label="Total Production", linewidth=2)
        plt.plot(years, exceed_production_tons, label="Exceedance Area Production", linewidth=2)
        plt.title(f"{basin} - {crop} (1986-2015)")
        plt.ylabel("Production (tons/year)")
        plt.xlabel("Year")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        fig_path = os.path.join(PlotDir, f"{basin}_{crop}_boundary_exceedance_production_1986-2015.png")
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
    plt.bar(x, avg_exceed_list, width, color=[crop_colors[c] for c in crop_labels], alpha=1.0, label="Exceedance Area Production")

    plt.xticks(x, crop_labels)
    plt.ylabel("Average Annual Production (tons/year)")
    plt.title(f"{basin} - Average Annual Production (1986-2015)")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    fig_path_bar = os.path.join(PlotDir, f"{basin}_summary_barchart.png")
    plt.savefig(fig_path_bar, dpi=300)
    plt.close()
    print(f"    ✅ Saved overlayed summary bar chart: {fig_path_bar}")


print("=== All processing complete ===")