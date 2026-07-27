import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# --------------------------------------------------
# Scenario paths
# --------------------------------------------------
scenarios = {
    "Current": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline",
    "Water": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_2_Sus_Irrigation",
    "Water + N": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_N_50",
    "Water + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_P_50",
    "Water + N + P": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Respect_NP_50"
}

# --------------------------------------------------
# Paths
# --------------------------------------------------
mask_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"

boundary_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Statistics/1_Boundary"

out_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V5_Plots/SM_Fig/Soilpool/IrriRainfed_Combined"
os.makedirs(out_dir, exist_ok=True)

studyareas = ["Rhine", "LaPlata", "Indus", "Yangtze"]
crops = ["winterwheat", "mainrice", "maize", "soybean"]

colors = ["#7A7878", "#458FBA", "#C06497", "#C89119", "#06956F"]

# --------------------------------------------------
# MAIN LOOP
# --------------------------------------------------
for basin in studyareas:
    for crop in crops:

        summary_file = os.path.join(
            "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline",
            f"{basin}_{crop}_summary.nc"
        )

        if not os.path.exists(summary_file):
            continue

        print(f"\nProcessing {basin} - {crop}")

        # --------------------------------------------------
        # Load spatial data
        # --------------------------------------------------
        ds = xr.open_dataset(summary_file)

        total_ha = ds["Total_HA"]
        bd = ds.get("bulk_density", None)

        ha_irri = ds["Irrigated_HA"]
        ha_rf = ds["Rainfed_HA"]

        valid = total_ha > 2500

        ha_irri = ha_irri.where(valid)
        ha_rf = ha_rf.where(valid)

        if bd is not None:
            bd = bd.where(valid)

        # --------------------------------------------------
        # Build mask dataframe
        # --------------------------------------------------
        if bd is not None:
            mask_df = xr.merge([
                ha_irri.rename("HA_irri"),
                ha_rf.rename("HA_rf"),
                bd
            ]).to_dataframe().reset_index()
        else:
            mask_df = xr.merge([
                ha_irri.rename("HA_irri"),
                ha_rf.rename("HA_rf")
            ]).to_dataframe().reset_index()

        mask_df = mask_df.dropna(subset=["HA_irri", "HA_rf"], how="all")

        # --------------------------------------------------
        # Load boundary threshold (NEW)
        # --------------------------------------------------
        boundary_file = os.path.join(
            boundary_dir,
            f"{basin}_crop-specific_boundaries.csv"
        )

        boundary_threshold = None

        if os.path.exists(boundary_file):

            bdf = pd.read_csv(boundary_file)

            # --------------------------------------------------
            # first column contains names (no header)
            # --------------------------------------------------
            name_col = bdf.columns[0]

            # match crop name in index-like column
            row = bdf[bdf[name_col].astype(str).str.strip() == crop]

            # fallback: sometimes naming mismatch
            if row.empty:
                row = bdf[bdf[name_col].astype(str).str.lower() == crop.lower()]

            # --------------------------------------------------
            # extract P boundary
            # --------------------------------------------------
            if not row.empty and "P [ktons]" in bdf.columns:

                p_val = row["P [ktons]"].values[0]

                if pd.notna(p_val):
                    boundary_threshold = float(p_val) * 0.25

        # --------------------------------------------------
        # storage
        # --------------------------------------------------
        pool_results = {s: [] for s in scenarios}
        runoff_results = {s: [] for s in scenarios}

        # --------------------------------------------------
        # scenario loop
        # --------------------------------------------------
        for sce, csv_dir in scenarios.items():

            csv_file = os.path.join(csv_dir, f"{basin}_{crop}_annual.csv")

            if not os.path.exists(csv_file):
                continue

            df = pd.read_csv(csv_file)

            for col in ["LabileP", "StableP", "P_surf", "P_sub"]:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors="coerce")

            df = df[(df["Year"] >= 2009) & (df["Year"] <= 2020)]

            merged = df.merge(
                mask_df,
                left_on=["Lat", "Lon"],
                right_on=["lat", "lon"],
                how="inner"
            )

            if merged.empty:
                continue

            for year, g in merged.groupby("Year"):

                # ==================================================
                # 1. SOIL P POOL (mmol/kg)
                # ==================================================
                if bd is not None:

                    w_irri = g["HA_irri"] * g["bulk_density"]
                    w_rf = g["HA_rf"] * g["bulk_density"]

                    num = 0
                    den = 0

                    if w_irri.sum() > 0:
                        num += ((g["LabileP"] + g["StableP"]) * w_irri).sum()
                        den += w_irri.sum()

                    if w_rf.sum() > 0:
                        num += ((g["LabileP"] + g["StableP"]) * w_rf).sum()
                        den += w_rf.sum()

                    pool_val = 30.97 * num / den if den > 0 else None

                else:

                    w_irri = g["HA_irri"]
                    w_rf = g["HA_rf"]

                    num = 0
                    den = 0

                    if w_irri.sum() > 0:
                        num += ((g["LabileP"] + g["StableP"]) * w_irri).sum()
                        den += w_irri.sum()

                    if w_rf.sum() > 0:
                        num += ((g["LabileP"] + g["StableP"]) * w_rf).sum()
                        den += w_rf.sum()

                    pool_val = 30.97 * num / den if den > 0 else None

                if pool_val is not None:
                    pool_results[sce].append((year, pool_val))

                # ==================================================
                # 2. P RUNOFF (kton)
                # ==================================================
                runoff_irri = ((g["P_surf"] + g["P_sub"]) * g["HA_irri"]).sum() * 1e-6
                runoff_rf = ((g["P_surf"] + g["P_sub"]) * g["HA_rf"]).sum() * 1e-6

                runoff_total = runoff_irri + runoff_rf

                runoff_results[sce].append((year, runoff_total))

        # --------------------------------------------------
        # PLOT
        # --------------------------------------------------
        fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

        # ----------------------------
        # Soil P pool
        # ----------------------------
        ax = axes[0]

        for i, (sce, vals) in enumerate(pool_results.items()):
            if len(vals) == 0:
                continue

            vals = pd.DataFrame(vals, columns=["Year", "Value"]).sort_values("Year")

            ax.plot(vals["Year"]+1, vals["Value"],label=sce,color=colors[i % len(colors)],lw=3)

        ax.grid(True, linestyle="--", alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_ylabel("Soil P pool [mg P/kg soil]", fontsize=15)

        # ----------------------------
        # P RUNOFF
        # ----------------------------
        ax = axes[1]

        # 🔴 NEW: boundary line
        if boundary_threshold is not None:
            ax.axhline(
                y=boundary_threshold,
                color="red",
                linestyle="--",
                linewidth=3,
                label="Safe boundary"
            )

        for i, (sce, vals) in enumerate(runoff_results.items()):
            if len(vals) == 0:
                continue

            vals = pd.DataFrame(vals, columns=["Year", "Value"]).sort_values("Year")

            ax.plot(vals["Year"]+1, vals["Value"],label=sce,color=colors[i % len(colors)],lw=2)

        # ax.set_title(f"{basin} - {crop} P runoff")
        # ax.set_xlabel("Year")
        ax.set_ylabel("P Runoff [kton P/yr]", fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.grid(True, linestyle="--", alpha=0.5)

        # axes[1].legend(loc="upper left", bbox_to_anchor=(1.05, 1))

        plt.tight_layout()

        out_file = os.path.join(
            out_dir,
            f"{basin}_{crop}_Ppool_runoff_combined.png"
        )

        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved: {out_file}")