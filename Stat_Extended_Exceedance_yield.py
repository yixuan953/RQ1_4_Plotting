# This script is used to calculate how many crop prouduction [Mt] is contributed by boundaries' exceedance cropland (irrigation, N runoff, P runoff) for each crop (and all major crops) and each basin

import pandas as pd
import numpy as np
import xarray as xr
import os

model_summary_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/Extended_Stat/4_Boundaries_Exceeded_Yield/2_Production"

data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

os.makedirs(output_dir, exist_ok=True)

Studyarea = ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]


def GetCropWNP(file_name, low_runoff_path, PA_path):

    with xr.open_dataset(file_name) as ds:
        mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)

    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]

    mask = mask.where(low_runoff != 1, np.nan)

    with xr.open_dataset(PA_path) as ds_pa:
        Total_PA = ds_pa["Physical_Area"] * mask
        Irri_PA = ds_pa["Physical_Area"] * ds_pa["Irrigated_Proportion"] * mask

    # --- Core variables ---
    Total_HA = ds["Total_HA"] * mask
    Irri_HA = ds["Irrigated_HA"] * mask

    T_Irr = ds["Total_irrigation_amount"] * mask
    S_Irr = ds["Sus_irrigation_amount"] * mask

    N = ds["N_Runoff"] * mask
    CN = ds["Crit_N_Runoff"] * mask

    P = ds["P_Runoff"] * mask
    CP = ds["Crit_P_Runoff"] * mask

    # --- NEW: yields ---
    Yield_Irr = ds["Avg_Yield_Irrigated"] * mask   # kg/ha
    Yield_RF  = ds["Avg_Yield_Rainfed"] * mask     # kg/ha

    return (Total_HA, Irri_HA,
            T_Irr, S_Irr,
            N, P, CN, CP,
            Total_PA, Irri_PA,
            Yield_Irr, Yield_RF)


# ================= MAIN LOOP =================

for basin in Studyarea:

    basin_results_prod = []
    basin_results_prod_pct = []

    basin_total_prod = 0.0
    basin_prod_sums = {
        "Irrigation": 0.0,
        "N": 0.0,
        "P": 0.0,
        "Irrigation & N": 0.0,
        "Irrigation & P": 0.0,
        "N & P": 0.0,
        "All three": 0.0
    }

    for crop in Croptypes:

        file_name = os.path.join(model_summary_dir, f"{basin}_{crop}_summary.nc")
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")

        if crop == "winterwheat":
            cropname = "WHEA"
        elif crop == "maize":
            cropname = "MAIZ"
        elif crop == "mainrice":
            cropname = "RICE"
        elif crop == "secondrice":
            cropname = "RICE"
        elif crop == "soybean":
            cropname = "SOYB"

        PA_path = os.path.join(
            data_dir,
            "2_StudyArea",
            basin,
            f"PhysicalArea/{cropname}_PA_05d.nc"
        )

        if os.path.exists(file_name) and os.path.exists(low_runoff_path) and os.path.exists(PA_path):

            (HA, Irri_HA,
             T_Irr, S_Irr,
             N, P, CN, CP,
             Total_PA, Irri_PA,
             Yield_Irr, Yield_RF) = GetCropWNP(file_name, low_runoff_path, PA_path)

            # --- masks ---
            w_exc = (T_Irr > S_Irr)
            n_exc = (N > CN)
            p_exc = (P > CP)

            Rainfed_HA = HA - Irri_HA

            # --- production (kg) ---
            Prod_Irr = Yield_Irr * Irri_HA
            Prod_RF = Yield_RF * Rainfed_HA
            Prod = Prod_Irr + Prod_RF

            total_prod = float(Prod.sum())

            if total_prod > 0:

                prod_exceeded = {
                    "Irrigation": float(Prod.where(w_exc).sum()),
                    "N": float(Prod.where(n_exc).sum()),
                    "P": float(Prod.where(p_exc).sum()),
                    "Irrigation & N": float(Prod.where(w_exc & n_exc).sum()),
                    "Irrigation & P": float(Prod.where(w_exc & p_exc).sum()),
                    "N & P": float(Prod.where(n_exc & p_exc).sum()),
                    "All three": float(Prod.where(w_exc & n_exc & p_exc).sum()),
                    "At least one": float(Prod.where(w_exc | n_exc | p_exc).sum()),
                    "No exceedance": float(Prod.where(~w_exc & ~n_exc & ~p_exc).sum()),
                }

                # --- accumulate basin totals ---
                basin_total_prod += total_prod
                for key in basin_prod_sums:
                    basin_prod_sums[key] += prod_exceeded[key]

                # --- store crop-level (Mt + %) ---
                row_mt = {"Crop type": crop}
                row_pct = {"Crop type": crop}

                for key in prod_exceeded:
                    row_mt[key] = prod_exceeded[key] / 1e9  # Mt
                    row_pct[key] = 100 * prod_exceeded[key] / total_prod

                basin_results_prod.append(row_mt)
                basin_results_prod_pct.append(row_pct)

    # --- All crops row ---
    if basin_total_prod > 0:

        all_mt = {"Crop type": "All Crops"}
        all_pct = {"Crop type": "All Crops"}

        for key in basin_prod_sums:
            all_mt[key] = basin_prod_sums[key] / 1e9
            all_pct[key] = 100 * basin_prod_sums[key] / basin_total_prod

        basin_results_prod.append(all_mt)
        basin_results_prod_pct.append(all_pct)

    # --- SAVE ---
    if basin_results_prod:

        df_mt = pd.DataFrame(basin_results_prod)
        df_pct = pd.DataFrame(basin_results_prod_pct)

        out_mt = os.path.join(output_dir, f"{basin}_boundary_exceedance_production_Mt.csv")
        out_pct = os.path.join(output_dir, f"{basin}_boundary_exceedance_production_pct.csv")

        df_mt.to_csv(out_mt, index=False)
        df_pct.to_csv(out_pct, index=False)

        print(f"Saved: {out_mt}")
        print(f"Saved: {out_pct}")