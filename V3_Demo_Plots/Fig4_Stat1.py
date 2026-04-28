# This script is used to summarize the crop production and N, P runoff changes under the 4 scenarios
# And save the output in the .csv file
# Column 1 - Crop production [ktons]
# Column 2 - N runoff [ktons]
# Column 3 - P runoff [ktons]

# Row 1 - Crop type: Wheat, Maize, Rice or Soybean
# Row 2 - Scenario: Baseline, Sus_Irri, Sus_Irri+Red_Fert, Sus_Irri+Red_Inc_Fert

import os
import xarray as xr
import numpy as np
import pandas as pd

# --- 1. Configuration ---
STUDY_AREAS = ["LaPlata", "Rhine", "Indus", "Yangtze"]
CROP_MAP = {
    "Wheat": ["winterwheat"],
    "Maize": ["maize"],
    "Rice": ["mainrice", "secondrice"], 
    "Soybean": ["soybean"]
}

SCENARIOS = {
    "Baseline": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline",
    "Sus_Irri": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/2_Sus_irrigation",
    "Sus_Irri+Red_Fert": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert",
    "Sus_Irri+Red_Inc_Fert": "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/4_Inc_fert"
}

DATA_DIR = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"
OUT_DIR  = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V3_Statistics/4_Red_Prod_Runoff"
os.makedirs(OUT_DIR, exist_ok=True)

# --- 2. Calculation Function ---

def get_total_production_runoff(folder, basin, crop_list):
    """Summarizes spatial data into basin-wide totals in ktons."""
    t_prod_sum = 0.0
    t_N_sum = 0.0
    t_P_sum = 0.0
    
    lr_path = os.path.join(DATA_DIR, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(lr_path) as ds_lr:
        runoff_mask = xr.where(ds_lr["Low_Runoff"].isnull(), 1, np.nan)

    for crop in crop_list:
        path = os.path.join(folder, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(path): 
            continue
            
        with xr.open_dataset(path) as ds:
            # Create mask: Only pixels with > 2500 HA and valid runoff
            crop_mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, np.nan)
            mask = crop_mask * runoff_mask

            # Production (assuming output is in tons, converting to ktons)
            prod = (ds["Avg_Yield_Irrigated"].fillna(0) * ds["Irrigated_HA"].fillna(0) + 
                    ds["Avg_Yield_Rainfed"].fillna(0) * ds["Rainfed_HA"].fillna(0)) * mask
            
            # Runoff (assuming output is in kg or tons, sum it up)
            # Note: Ensure N_Runoff and P_Runoff in your NC files are MASS (kg or tons), not concentrations
            n_run = ds["N_Runoff"] * mask
            p_run = ds["P_Runoff"] * mask

            t_prod_sum += prod.sum().values
            t_N_sum += n_run.sum().values
            t_P_sum += p_run.sum().values
            
    # Convert to ktons 
    return t_prod_sum / 1e6, t_N_sum / 1e6, t_P_sum / 1e6

# --- 3. Execution and CSV Export ---

all_data = []

for basin in STUDY_AREAS:
    for crop_name, crop_files in CROP_MAP.items():
        for scen_name, scen_path in SCENARIOS.items():
            
            prod_kton, n_kton, p_kton = get_total_production_runoff(scen_path, basin, crop_files)
            
            # Append result row
            all_data.append({
                "Basin": basin,
                "Crop": crop_name,
                "Scenario": scen_name,
                "Production_ktons": prod_kton,
                "N_Runoff_ktons": n_kton,
                "P_Runoff_ktons": p_kton
            })

# Create DataFrame
df = pd.DataFrame(all_data)

# Save to CSV
csv_path = os.path.join(OUT_DIR, "Production_Runoff_Summary.csv")
df.to_csv(csv_path, index=False)

print(f"Summary successfully saved to: {csv_path}")