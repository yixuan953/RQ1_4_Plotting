# This code is used to summarice if the regional boundary is met or not

import os
import xarray as xr
import numpy as np

Studyarea =  ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes =  ["winterwheat", "maize", "mainrice", "secondrice", "soybean"] 

for basin in Studyarea:
    output_file = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert/{basin}_all_crop_summary.nc"
    range_nc = os.path.join("/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea", basin, "range.nc")
    
    with xr.open_dataset(range_nc) as ds_range:
        template = ds_range["mask"]
        
    # Initialize totals as zeros with the EXACT same coordinates as the template
    total_n = xr.full_like(template, 0.0, dtype=np.float32)
    total_p = xr.full_like(template, 0.0, dtype=np.float32)
    total_crit_n = xr.full_like(template, 0.0, dtype=np.float32)
    total_crit_p = xr.full_like(template, 0.0, dtype=np.float32)

    for crop in Croptypes:
        input_file = f"/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/3_Red_fert/{basin}_{crop}_summary.nc"

        if not os.path.exists(input_file):
            print(f"Skipping: {crop} in {basin} (File not found)")
            continue

        with xr.open_dataset(input_file) as ds:
            # Reindex ensures the crop data matches the template coordinates exactly
            # This solves the 0.0000001 difference errors
            ds_aligned = ds.reindex_like(template, method="nearest")
            
            mask = ds_aligned['Basin_mask'].where(ds_aligned['Total_HA'] > 2500)

            # Use += with fillna to accumulate values
            total_n += (ds_aligned['N_Runoff'] * mask).fillna(0)
            total_p += (ds_aligned['P_Runoff'] * mask).fillna(0)
            total_crit_n += (ds_aligned['Crit_N_Runoff'] * mask).fillna(0)
            total_crit_p += (ds_aligned['Crit_P_Runoff'] * mask).fillna(0)

    # Convert zeros back to NaN where no data was added (using the template's own mask)
    # This keeps the basin boundaries clean
    final_mask = template != 0 
    
    ds_out = xr.Dataset({
        "N_Runoff": total_n.where(final_mask, np.nan),
        "P_Runoff": total_p.where(final_mask, np.nan),
        "Crit_N_Runoff": total_crit_n.where(final_mask, np.nan),
        "Crit_P_Runoff": total_crit_p.where(final_mask, np.nan),
    })

    ds_out.to_netcdf(output_file)
    print(f"Successfully saved: {output_file}")