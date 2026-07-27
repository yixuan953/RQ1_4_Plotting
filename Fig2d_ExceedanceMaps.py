import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# --------------------------------------------------
# MODE: "Irrigated" or "Rainfed"
# --------------------------------------------------
mode = "Irrigated"

# --------------------------------------------------
# Scenario paths
# --------------------------------------------------
scenarios_irrigated = {
    "Current": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline",
    "Water": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_2_Sus_Irrigation",
    "Water + N": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_N_50",
    "Water + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_P_50",
    "Water + N + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_NP_50"
}

scenarios_rainfed = {
    "Current": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline_rainfed",
    "Water": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline_rainfed",
    "Water + N": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed/Respect_N_50",
    "Water + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed/Respect_P_50",
    "Water + N + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed/Respect_NP_50"
}

scenarios = scenarios_irrigated if mode == "Irrigated" else scenarios_rainfed

# --------------------------------------------------
# Paths
# --------------------------------------------------
mask_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"

out_dir = f"/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/Extended_Fig/Fig7/{mode}"
os.makedirs(out_dir, exist_ok=True)

studyareas = ["Rhine"] # ["LaPlata", "Yangtze", "Indus", "Rhine"]
crops = ["winterwheat"] # ["mainrice", "secondrice", "winterwheat", "soybean", "maize"]

# --------------------------------------------------
# Helper
# --------------------------------------------------
def fix_crop(crop):
    return "winterwheat" if crop == "wheat" else crop

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

# --------------------------------------------------
# MAIN LOOP
# --------------------------------------------------
for basin in studyareas:
    for crop in crops:

        crop_name = fix_crop(crop)

        summary_file = os.path.join(
            "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline",
            f"{basin}_{crop_name}_summary.nc"
        )

        if not os.path.exists(summary_file):
            continue

        print(f"\nProcessing {basin} - {crop}")

        # --------------------------------------------------
        # Load spatial summary (HA + BD + masks)
        # --------------------------------------------------
        ds = xr.open_dataset(summary_file)

        total_ha = ds["Total_HA"]
        bd = ds.get("bulk_density", None)

        if mode == "Irrigated":
            ha = ds["Irrigated_HA"]
        else:
            ha = ds["Rainfed_HA"]

        # apply filter
        valid = total_ha > 2500

        ha = ha.where(valid)

        if bd is not None:
            bd = bd.where(valid)

        # dataframe for merging
        if bd is not None:
            mask_df = xr.merge([ha, bd]).to_dataframe().reset_index()
        else:
            mask_df = ha.to_dataframe(name="HA").reset_index()

        mask_df = mask_df.dropna(subset=["HA"])

        # --------------------------------------------------
        # storage
        # --------------------------------------------------
        pool_results = {s: [] for s in scenarios}
        loss_results = {s: [] for s in scenarios}

        # --------------------------------------------------
        # scenario loop
        # --------------------------------------------------
        for sce, csv_dir in scenarios.items():

            csv_file = os.path.join(csv_dir, f"{basin}_{crop}_annual.csv")

            if not os.path.exists(csv_file):
                continue

            df = pd.read_csv(csv_file)

            for col in ["LabileP", "StableP", "P_surf", "Psub"]:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors="coerce")

            df = df[(df["Year"] >= 2005) & (df["Year"] <= 2019)]

            merged = df.merge(
                mask_df,
                left_on=["Lat", "Lon"],
                right_on=["lat", "lon"],
                how="inner"
            )

            if merged.empty:
                continue

            for year, g in merged.groupby("Year"):

                # ----------------------------
                # Soil P pools (HA × BD)
                # ----------------------------
                if bd is not None:
                    w_pool = g["HA"] * g["bulk_density"]
                else:
                    w_pool = g["HA"]

                denom_pool = w_pool.sum()
                if denom_pool > 0:
                    pool_val = (
                        (g["LabileP"] + g["StableP"]) * w_pool
                    ).sum() / denom_pool

                    pool_results[sce].append((year, pool_val))

                # ----------------------------
                # P losses (HA only)
                # ----------------------------
                w_loss = g["HA"]
                denom_loss = w_loss.sum()

                if denom_loss > 0:
                    loss_val = (
                        (g["P_surf"] + g["Psub"]) * w_loss
                    ).sum() / denom_loss

                    loss_results[sce].append((year, loss_val))

        # --------------------------------------------------
        # PLOT
        # --------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # ----------------------------
        # Subplot 1: Soil P pools
        # ----------------------------
        ax = axes[0]

        for i, (sce, vals) in enumerate(pool_results.items()):
            if len(vals) == 0:
                continue

            vals = pd.DataFrame(vals, columns=["Year", "Value"]).sort_values("Year")

            ax.plot(vals["Year"], vals["Value"],
                    label=sce,
                    color=colors[i % len(colors)],
                    lw=2)

        ax.set_title(f"{basin} - {crop} Soil P pools")
        ax.set_xlabel("Year")
        ax.set_ylabel("Labile + Stable P")

        # ----------------------------
        # Subplot 2: P losses
        # ----------------------------
        ax = axes[1]

        for i, (sce, vals) in enumerate(loss_results.items()):
            if len(vals) == 0:
                continue

            vals = pd.DataFrame(vals, columns=["Year", "Value"]).sort_values("Year")

            ax.plot(vals["Year"], vals["Value"],
                    label=sce,
                    color=colors[i % len(colors)],
                    lw=2)

        ax.set_title(f"{basin} - {crop} P losses")
        ax.set_xlabel("Year")
        ax.set_ylabel("P_surf + Psub")

        axes[1].legend(loc="upper left", bbox_to_anchor=(1.05, 1))

        plt.tight_layout()

        out_file = os.path.join(out_dir, f"{basin}_{crop}_{mode}_Ppool_loss.png")
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved: {out_file}")