import os
import xarray as xr
import numpy as np

# Configuration
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
InputCrops = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

# Paths
org_fert_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
final_fert_dir_irrigated = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/4_Fertilization_Red/4_1_Reduced_Fert/Irrigated/Inc_prop"
final_fert_dir_rainfed = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/4_Fertilization_Red/4_1_Reduced_Fert/Rainfed/Inc_prop"
summary_file_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/4_Inc_fert"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/1_Fertilizer_input"

def calculate_total_fertilizer(fert_file, HA, cropname, applied_mask):
    """Calculates total N and P in kg, applying the spatial mask and HA filter."""
    with xr.open_dataset(fert_file) as fert_ds:
        # Select 2015 data
        Urea_N = fert_ds[f"{cropname}_Urea_N_application_rate"].sel(year=2015)
        Inorg_N = fert_ds[f"{cropname}_Inorg_N_application_rate"].sel(year=2015)
        Manure_N = fert_ds[f"{cropname}_Manure_N_application_rate"].sel(year=2015)
        Inorg_P = fert_ds[f"{cropname}_P_application_rate"].sel(year=2015)
        Manure_P = fert_ds[f"{cropname}_Manure_P_application_rate"].sel(year=2015)
        
        # Calculate total inputs (kg) applied only to valid mask areas
        total_N_kg = (Urea_N + Inorg_N + Manure_N) * HA * applied_mask
        total_P_kg = (Inorg_P + Manure_P) * HA * applied_mask
        
    return total_N_kg, total_P_kg 

for basin in Studyareas:
    print(f"\n{'='*60}\nProcessing Basin: {basin}\n{'='*60}")
    
    # Load basin mask
    basin_mask_file = os.path.join(org_fert_dir, basin, "range.nc")
    if not os.path.exists(basin_mask_file):
        print(f"Mask missing for {basin}, skipping...")
        continue
    
    ds_basin_mask = xr.open_dataset(basin_mask_file)
    mask = ds_basin_mask["mask"].where(ds_basin_mask["mask"] == 1.0)

    # Store rice datasets for Yangtze summation
    yangtze_rice_list = []

    for crop in InputCrops:
        summary_file = os.path.join(summary_file_dir, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(summary_file):
            continue

        ds_summary = xr.open_dataset(summary_file)
        Total_HA = ds_summary["Total_HA"]
        Irrigated_HA = ds_summary["Irrigated_HA"]
        Rainfed_HA = ds_summary["Rainfed_HA"]
        
        # Final mask: Basin boundary + cells where Harvested Area > 2500 ha
        Basin_mask = mask.where(Total_HA > 2500, np.nan)

        # Map internal crop names to file naming conventions
        mapping = {
            "mainrice": ("Rice", "Rice"),
            "secondrice": ("Rice", "Secondrice"),
            "soybean": ("Soybean", "Soybean"),
            "maize": ("Maize", "Maize"),
            "winterwheat": ("Wheat", "Wheat")
        }
        cropname, cropname2 = mapping[crop]

        # File Paths
        org_nc = os.path.join(org_fert_dir, f"{basin}/Fertilization/{basin}_{cropname}_Fert_2005-2020_FixRate.nc")
        irrig_nc = os.path.join(final_fert_dir_irrigated, f"{basin}_{cropname2}_Fert_2005-2020_FixRate.nc")
        rain_nc = os.path.join(final_fert_dir_rainfed, f"{basin}_{cropname2}_Fert_2005-2020_FixRate.nc")

        if not all(os.path.exists(f) for f in [org_nc, irrig_nc, rain_nc]):
            print(f"Skipping {crop}: missing input files.")
            continue

        # Core Calculations
        org_N, org_P = calculate_total_fertilizer(org_nc, Total_HA, cropname, Basin_mask)
        rain_N, rain_P = calculate_total_fertilizer(rain_nc, Rainfed_HA, cropname, Basin_mask)
        irrig_N, irrig_P = calculate_total_fertilizer(irrig_nc, Irrigated_HA, cropname, Basin_mask)

        final_N = (rain_N.fillna(0) + irrig_N.fillna(0)) * Basin_mask
        final_P = (rain_P.fillna(0) + irrig_P.fillna(0)) * Basin_mask

        # Unit Conversion (kg to ktons)
        conv = 0.000001
        b_org_N, b_final_N = org_N.sum().item() * conv, final_N.sum().item() * conv
        b_org_P, b_final_P = org_P.sum().item() * conv, final_P.sum().item() * conv

        print(f"Crop: {crop:<12} | N: {b_org_N:>8.2f} -> {b_final_N:>8.2f} ktons | P: {b_org_P:>7.2f} -> {b_final_P:>7.2f} ktons")

        # Create Individual Crop Dataset
        ds_output = xr.Dataset({
            "org_fert_total_N": org_N,
            "org_fert_total_P": org_P,
            "final_fert_total_N": final_N,
            "final_fert_total_P": final_P,
            "N_change": final_N - org_N,
            "N_change_kgperha": (final_N - org_N)/Total_HA,
            "P_change": final_P - org_P,
            "P_change_kgperha": (final_P - org_P)/Total_HA
        })

        # Save individual crop
        output_nc = os.path.join(output_dir, f"{basin}_{crop}_Fert_Change.nc")
        ds_output.to_netcdf(output_nc)

        # Yangtze Rice Aggregation
        if basin == "Yangtze" and crop in ["mainrice", "secondrice"]:
            yangtze_rice_list.append(ds_output)

    # After processing all crops in Yangtze, sum the rice
    if basin == "Yangtze" and len(yangtze_rice_list) == 2:
        print("-" * 60)
        total_rice_ds = yangtze_rice_list[0] + yangtze_rice_list[1]
        
        tr_org_N = total_rice_ds["org_fert_total_N"].sum().item() * conv
        tr_final_N = total_rice_ds["final_fert_total_N"].sum().item() * conv
        
        print(f"Combined Yangtze Rice | N: {tr_org_N:>8.2f} -> {tr_final_N:>8.2f} ktons")
        
        combined_output = os.path.join(output_dir, "Yangtze_totalrice_Fert_Change.nc")
        total_rice_ds.to_netcdf(combined_output)
        print(f"Saved Combined Rice: {combined_output}")