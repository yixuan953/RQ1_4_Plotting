# This code is used to plot the N, P uptakes and losses through water flow under the following conditions: 
# 1) Unsustainble irrigation
# 2) Sustainable irrigatin
# 3) Rainfed

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Input and output directions
DataDir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
OutputDir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs"

Basins = ["Indus"]  # ["Indus", "Rhine", "LaPlata", "Yangtze"] 
CropTypes = ["winterwheat"] # ["mainrice", "maize", "winterwheat", "soybean", "secondrice"]

PlotDir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Part2"

# Years for statistics 
start_year = 2005
end_year = 2015

# Pick colors for uptakes and losses
bar_colors = {
    "Crop uptake": "#E69F00",
    "Runoff": "#56B4E9",
    "Leaching": "#0072B2"
}

for basin in Basins:
    print(f"=== Processing basin: {basin} ===")

    for crop in CropTypes:
        print(f"> {crop}")

        # Read model output data
        output_nc_path_unsus_irrigated = os.path.join(OutputDir, f"Output_Unsus_Irrigation/{basin}_{crop}_annual.nc")
        output_nc_path_sus_irrigated = os.path.join(OutputDir, f"Output_Sus_Irrigation/{basin}_{crop}_annual.nc")
        output_nc_path_rainfed = os.path.join(OutputDir, f"Output_Rainfed/{basin}_{crop}_annual.nc")
        # Check path
        if not os.path.exists(output_nc_path_unsus_irrigated):
            print(f"Missing output file: {output_nc_path_unsus_irrigated}")
            continue
        if not os.path.exists(output_nc_path_sus_irrigated):
            print(f"Missing output file: {output_nc_path_sus_irrigated}")
            continue
        if not os.path.exists(output_nc_path_rainfed):
            print(f"Missing output file: {output_nc_path_rainfed}")
            continue
        output_nc_unsus_irrigated = xr.open_dataset(output_nc_path_unsus_irrigated).sel(time=slice(start_year, end_year))
        output_nc_sus_irrigated = xr.open_dataset(output_nc_path_sus_irrigated).sel(time=slice(start_year, end_year))
        output_nc_rainfed = xr.open_dataset(output_nc_path_rainfed).sel(time=slice(start_year, end_year))

        N_fert_input = output_nc_unsus_irrigated["N_fert"].fillna(0) + output_nc_unsus_irrigated["N_surf"].fillna(0) + output_nc_unsus_irrigated["NH3"].fillna(0) + output_nc_unsus_irrigated["NOx"].fillna(0) + output_nc_unsus_irrigated["N2O"].fillna(0)
        P_fert_input = output_nc_unsus_irrigated["P_fert"]
        
        # Compute N, P losses for irrigated cropland
        # Unsustainable
        N_uptake_unsus_irrigated = output_nc_unsus_irrigated["N_uptake"]
        N_leach_unsus_irrigated = output_nc_unsus_irrigated["N_leach"]
        N_runoff_unsus_irrigated = output_nc_unsus_irrigated["N_surf"].fillna(0) + output_nc_unsus_irrigated["N_sub"].fillna(0)
        P_uptake_unsus_irrigated = output_nc_unsus_irrigated["P_uptake"]
        P_leach_unsus_irrigated = output_nc_unsus_irrigated["P_leach"]
        P_runoff_unsus_irrigated = output_nc_unsus_irrigated["P_surf"].fillna(0) + output_nc_unsus_irrigated["P_sub"].fillna(0)
        N_decomp_unsus_irrigated = output_nc_unsus_irrigated["N_decomp"]
        P_decomp_unsus_irrigated = output_nc_unsus_irrigated["P_decomp"]
        # Sustainable 
        N_uptake_sus_irrigated = output_nc_sus_irrigated["N_uptake"]
        N_leach_sus_irrigated = output_nc_sus_irrigated["N_leach"]
        N_runoff_sus_irrigated = output_nc_sus_irrigated["N_surf"].fillna(0) + output_nc_sus_irrigated["N_sub"].fillna(0)
        P_uptake_sus_irrigated = output_nc_sus_irrigated["P_uptake"]
        P_leach_sus_irrigated = output_nc_sus_irrigated["P_leach"]
        P_runoff_sus_irrigated = output_nc_sus_irrigated["P_surf"].fillna(0) + output_nc_sus_irrigated["P_sub"].fillna(0)
        N_decomp_sus_irrigated = output_nc_sus_irrigated["N_decomp"]
        P_decomp_sus_irrigated = output_nc_sus_irrigated["P_decomp"]
        # Load harvested area for irrigated 
        HA_path_irrigated = os.path.join(DataDir, basin, "Mask", f"{crop}_HA_Irrigated.nc")
        if not os.path.exists(HA_path_irrigated):
            print(f"Missing mask file: {HA_path_irrigated}, skipping crop.")
            continue
        HA_nc_irrigated = xr.open_dataset(HA_path_irrigated)
        HA_irrigated = HA_nc_irrigated["area_total"] if "area_total" in HA_nc_irrigated else list(HA_nc_irrigated.data_vars.values())[0]

        # Compute N, P losses for rainfed cropland
        N_uptake_rainfed = output_nc_rainfed["N_uptake"]
        N_leach_rainfed = output_nc_rainfed["N_leach"]
        N_runoff_rainfed = output_nc_rainfed["N_surf"].fillna(0) + output_nc_rainfed["N_sub"].fillna(0)
        P_uptake_rainfed = output_nc_rainfed["P_uptake"]
        P_leach_rainfed = output_nc_rainfed["P_leach"]
        P_runoff_rainfed = output_nc_rainfed["P_surf"].fillna(0) + output_nc_rainfed["P_sub"].fillna(0)
        N_decomp_rainfed = output_nc_rainfed["N_decomp"]
        P_decomp_rainfed = output_nc_rainfed["P_decomp"]
        # Load harvested area for irrigated 
        HA_path_rainfed = os.path.join(DataDir, basin, "Mask", f"{crop}_HA_rainfed.nc")
        if not os.path.exists(HA_path_rainfed):
            print(f"Missing mask file: {HA_path_rainfed}, skipping crop.")
            continue
        HA_nc_rainfed = xr.open_dataset(HA_path_rainfed)
        HA_rainfed = HA_nc_rainfed["area_total"] if "area_total" in HA_nc_rainfed else list(HA_nc_rainfed.data_vars.values())[0]
        HA_total = HA_rainfed.fillna(0) + HA_irrigated.fillna(0)

        HA_total = xr.DataArray(HA_total, dims=("lat", "lon"))
        HA_total = HA_total.interp(lat=N_fert_input.lat, lon=N_fert_input.lon, method="nearest")
        valid_HA = xr.where(HA_total > 2500, 1, 0)

        # Get the weighted average N, P decomposition, fetilizer input, uptakes and losses for each basin
        N_fert_WA = (
            (N_fert_input * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_fert_WA = (
            (P_fert_input * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        # Unsustainable irrigated
        N_uptake_unsus_irri_WA = (
            (N_uptake_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_runoff_unsus_irri_WA = (
            (N_runoff_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_leach_unsus_irri_WA = (
            (N_leach_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_decomp_unsus_irri_WA = (
            (N_decomp_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )       
        P_uptake_unsus_irri_WA = (
            (P_uptake_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_runoff_unsus_irri_WA = (
            (P_runoff_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_leach_unsus_irri_WA = (
            (P_leach_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_decomp_unsus_irri_WA = (
            (P_decomp_unsus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )

        # Sustainable irrigated
        N_uptake_sus_irri_WA = (
            (N_uptake_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_runoff_sus_irri_WA = (
            (N_runoff_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_leach_sus_irri_WA = (
            (N_leach_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_decomp_sus_irri_WA = (
            (N_decomp_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )       
        P_uptake_sus_irri_WA = (
            (P_uptake_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_runoff_sus_irri_WA = (
            (P_runoff_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_leach_sus_irri_WA = (
            (P_leach_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_decomp_sus_irri_WA = (
            (P_decomp_sus_irrigated * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        ) 

        # Rainfed - CORRECTED VARIABLE NAMES
        N_uptake_rainfed_WA = (
            (N_uptake_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_runoff_rainfed_WA = (
            (N_runoff_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_leach_rainfed_WA = (
            (N_leach_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        N_decomp_rainfed_WA = (
            (N_decomp_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )       
        P_uptake_rainfed_WA = (
            (P_uptake_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_runoff_rainfed_WA = (
            (P_runoff_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_leach_rainfed_WA = (
            (P_leach_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )
        P_decomp_rainfed_WA = (
            (P_decomp_rainfed * valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
            / (valid_HA * HA_irrigated).sum(dim=["lat", "lon"], skipna=True)
        )      

        #  ===========  Start plotting ===========
        def to_np(x):
            return x.values if isinstance(x, xr.DataArray) else x
        
        # === N variables ===
        N_uptake_unsus = to_np(N_uptake_unsus_irri_WA)
        N_runoff_unsus = to_np(N_runoff_unsus_irri_WA)
        N_leach_unsus = to_np(N_leach_unsus_irri_WA)
        N_decomp_unsus = to_np(N_decomp_unsus_irri_WA)

        N_uptake_sus = to_np(N_uptake_sus_irri_WA)
        N_runoff_sus = to_np(N_runoff_sus_irri_WA)
        N_leach_sus = to_np(N_leach_sus_irri_WA)
        N_decomp_sus = to_np(N_decomp_sus_irri_WA)

        N_uptake_rain = to_np(N_uptake_rainfed_WA)
        N_runoff_rain = to_np(N_runoff_rainfed_WA)
        N_leach_rain = to_np(N_leach_rainfed_WA)
        N_decomp_rain = to_np(N_decomp_rainfed_WA)

        N_fert = to_np(N_fert_WA)

        # === P variables ===
        P_uptake_unsus = to_np(P_uptake_unsus_irri_WA)
        P_runoff_unsus = to_np(P_runoff_unsus_irri_WA)
        P_leach_unsus = to_np(P_leach_unsus_irri_WA)
        P_decomp_unsus = to_np(P_decomp_unsus_irri_WA)

        P_uptake_sus = to_np(P_uptake_sus_irri_WA)
        P_runoff_sus = to_np(P_runoff_sus_irri_WA)
        P_leach_sus = to_np(P_leach_sus_irri_WA)
        P_decomp_sus = to_np(P_decomp_sus_irri_WA)

        P_uptake_rain = to_np(P_uptake_rainfed_WA)
        P_runoff_rain = to_np(P_runoff_rainfed_WA)
        P_leach_rain = to_np(P_leach_rainfed_WA)
        P_decomp_rain = to_np(P_decomp_rainfed_WA)

        P_fert = to_np(P_fert_WA)

        # ================ Define groups
        systems = ["Unsus. Irrig.", "Sus. Irrig.", "Rainfed"]
        x = np.arange(len(systems))
        width = 0.3  # bar width per group

        # --- Helper functions ---
        def extract_array(data):
            """Return a 1D numpy array from DataArray or Dataset."""
            if isinstance(data, xr.Dataset):
                data = list(data.data_vars.values())[0]
            if isinstance(data, xr.DataArray):
                return data.values
            else:
                # Already a numpy array
                return data

        def extract_mean(data):
            """Return mean as float from DataArray, Dataset, NumPy array, or list."""
            if isinstance(data, xr.Dataset):
                # Extract first variable if Dataset
                data = list(data.data_vars.values())[0]
            if isinstance(data, xr.DataArray):
                return float(data.mean(skipna=True).values)
            else:
                # NumPy array or list — just use np.nanmean to ignore NaNs
                return float(np.nanmean(data))

        # --- Create figure ---
        fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)

        # ======================================================
        #                  SUBPLOT 1: NITROGEN
        # ======================================================
        ax = axes[0]

        # --- Bar data ---
        for i, (uptake, runoff, leach) in enumerate([
            (N_uptake_unsus, N_runoff_unsus, N_leach_unsus),
            (N_uptake_sus, N_runoff_sus, N_leach_sus),
            (N_uptake_rain, N_runoff_rain, N_leach_rain),
        ]):
            uptake_mean = extract_mean(uptake)
            runoff_mean = extract_mean(runoff)
            leach_mean  = extract_mean(leach)

            ax.bar(x[i], leach_mean, color=bar_colors["Leaching"], width=width)
            ax.bar(x[i], runoff_mean, bottom=leach_mean, color=bar_colors["Runoff"], width=width)
            ax.bar(x[i], uptake_mean, bottom=leach_mean + runoff_mean, color=bar_colors["Crop uptake"], width=width)

        # --- Boxplot for total N fluxes ---
        total_N = np.vstack([
            extract_array(N_uptake_unsus) + extract_array(N_runoff_unsus) + extract_array(N_leach_unsus),
            extract_array(N_uptake_sus) + extract_array(N_runoff_sus) + extract_array(N_leach_sus),
            extract_array(N_uptake_rain) + extract_array(N_runoff_rain) + extract_array(N_leach_rain)
        ]).T

        ax.boxplot(total_N, positions=x, widths=0.08, patch_artist=True,
                boxprops=dict(facecolor='white', color='black'),
                medianprops=dict(color='black'))

        # --- Fertilizer & decomposition lines ---
        ax.axhline(extract_mean(N_fert), color='red', linestyle='--', label='Fertilizer input')
        ax.plot(x, [extract_mean(N_decomp_unsus), extract_mean(N_decomp_sus), extract_mean(N_decomp_rain)],
                color='gray', linestyle='--', marker='o', label='Decomposition')

        # --- Axis & labels ---
        ax.set_xticks(x)
        ax.set_xticklabels(systems)
        ax.set_ylim(0, 150)
        ax.set_ylabel("kg N ha⁻¹ yr⁻¹")
        ax.set_title("N")
        
        # --- Create custom legend handles for N subplot ---
        # legend_handles_N = [
        #     mpatches.Patch(color=bar_colors["Crop uptake"], label="Crop uptake"),
        #     mpatches.Patch(color=bar_colors["Runoff"], label="Runoff"),
        #     mpatches.Patch(color=bar_colors["Leaching"], label="Leaching"),
        #     plt.Line2D([0], [0], color='red', linestyle='--', label='Fertilizer input'),
        #     plt.Line2D([0], [0], color='gray', linestyle='--', marker='o', label='Decomposition')
        # ]
        # ax.legend(handles=legend_handles_N, loc='upper right', frameon=True)

        # ======================================================
        #                  SUBPLOT 2: PHOSPHORUS
        # ======================================================
        ax = axes[1]

        for i, (uptake, runoff, leach) in enumerate([
            (P_uptake_unsus, P_runoff_unsus, P_leach_unsus),
            (P_uptake_sus, P_runoff_sus, P_leach_sus),
            (P_uptake_rain, P_runoff_rain, P_leach_rain),
        ]):
            uptake_mean = extract_mean(uptake)
            runoff_mean = extract_mean(runoff)
            leach_mean  = extract_mean(leach)

            ax.bar(x[i], leach_mean, color=bar_colors["Leaching"], width=width)
            ax.bar(x[i], runoff_mean, bottom=leach_mean, color=bar_colors["Runoff"], width=width)
            ax.bar(x[i], uptake_mean, bottom=leach_mean + runoff_mean, color=bar_colors["Crop uptake"], width=width)

        # --- Boxplot for total P fluxes ---
        total_P = np.vstack([
            extract_array(P_uptake_unsus) + extract_array(P_runoff_unsus) + extract_array(P_leach_unsus),
            extract_array(P_uptake_sus) + extract_array(P_runoff_sus) + extract_array(P_leach_sus),
            extract_array(P_uptake_rain) + extract_array(P_runoff_rain) + extract_array(P_leach_rain)
        ]).T

        ax.boxplot(total_P, positions=x, widths=0.08, patch_artist=True,
                boxprops=dict(facecolor='white', color='black'),
                medianprops=dict(color='black'))

        # --- Fertilizer & decomposition lines ---
        ax.axhline(extract_mean(P_fert), color='red', linestyle='--', label='Fertilizer input')
        ax.plot(x, [extract_mean(P_decomp_unsus), extract_mean(P_decomp_sus), extract_mean(P_decomp_rain)],
                color='gray', linestyle='--', marker='o', label='Decomposition')

        # --- Axis & labels ---
        ax.set_xticks(x)
        ax.set_xticklabels(systems)
        ax.set_ylim(0, 40)
        ax.set_ylabel("kg P ha⁻¹ yr⁻¹")
        ax.set_title("P")
        
        # --- Create custom legend handles for P subplot ---
        legend_handles_P = [
            mpatches.Patch(color=bar_colors["Crop uptake"], label="Crop uptake"),
            mpatches.Patch(color=bar_colors["Runoff"], label="Runoff"),
            mpatches.Patch(color=bar_colors["Leaching"], label="Leaching"),
            plt.Line2D([0], [0], color='red', linestyle='--', label='Fertilizer input'),
            plt.Line2D([0], [0], color='gray', linestyle='--', marker='o', label='Decomposition')
        ]
        ax.legend(handles=legend_handles_P, loc='upper right', frameon=True)

        plt.tight_layout()
        fig_path_bar = os.path.join(PlotDir, f"{basin}_{crop}_Fig2.png")
        plt.savefig(fig_path_bar, dpi=300)
        plt.close()
        print(f"✅ Saved bar charts: {fig_path_bar}")