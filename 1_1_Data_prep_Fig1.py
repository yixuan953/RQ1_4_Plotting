import os
import xarray as xr
import numpy as np

start_year = 2010
end_year = 2019

Studyarea =  ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes =  ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary"

# Input directory 1 - Simulated yield and losses
# Raifed_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline_rainfed"
Irrigated_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_2_Sus_Irrigation"  # Sustainable irrigation
# Irrigated_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_1_Baseline"

# Raifed_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Rainfed/Red_org"
# Irrigated_baseline_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/3_Scenarios/2_3_Sus_Irri_Red_Fert/Red_org"

# Input directory 2 - Model input data, including, Harvested area from SPAM2010, Irrigation from VIC-WUR
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"

# Input directory 3 - Calculated critical N, P runoff
crit_loss_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/2_Critical_NP_losses/Method3"

for basin in Studyarea:
    
    # Load basin range
    basin_mask_file = os.path.join(data_dir, basin, f"range.nc")
    ds_basin_mask = xr.open_dataset(basin_mask_file)
    mask = ds_basin_mask["mask"]
    mask = mask.where(mask == 1.0, np.nan)

    for crop in Croptypes:

        # Load simulated yield and runoff [kg/ha]
        Rainfed_baseline_file = os.path.join(Raifed_baseline_dir, f"{basin}_{crop}_annual.nc")
        Irrigated_baseline_file = os.path.join (Irrigated_baseline_dir, f"{basin}_{crop}_annual.nc")
        if not os.path.exists (Rainfed_baseline_file) or not os.path.exists (Irrigated_baseline_file):
            print (f"Missing {crop} for basin {basin}")
            continue

        ds_rainfed_sim = xr.open_dataset(Rainfed_baseline_file)
        N_Runoff_Rainfed_all = ds_rainfed_sim["N_Runoff"].sel(year=slice(start_year, end_year))
        P_Runoff_Rainfed_all = ds_rainfed_sim["P_Runoff"].sel(year=slice(start_year, end_year))
        Yield_Rainfed = ds_rainfed_sim["Yield"].sel(year=slice(start_year, end_year))
        N_Runoff_Rainfed =  N_Runoff_Rainfed_all.mean(dim = "year", skipna=True)
        P_Runoff_Rainfed =  P_Runoff_Rainfed_all.mean(dim = "year", skipna=True)
        Avg_Yield_Rainfed = Yield_Rainfed.mean(dim = "year", skipna=True)
        
        ds_Irrigated_sim = xr.open_dataset(Irrigated_baseline_file)
        N_Runoff_Irrigated_all = ds_Irrigated_sim["N_Runoff"].sel(year=slice(start_year, end_year))
        P_Runoff_Irrigated_all = ds_Irrigated_sim["P_Runoff"].sel(year=slice(start_year, end_year))
        Yield_Irrigated = ds_Irrigated_sim["Yield"].sel(year=slice(start_year, end_year))
        N_Runoff_Irrigated =  N_Runoff_Irrigated_all.mean(dim = "year", skipna=True)
        P_Runoff_Irrigated =  P_Runoff_Irrigated_all.mean(dim = "year", skipna=True)
        Avg_Yield_Irrigated = Yield_Irrigated.mean(dim = "year", skipna=True)

        # Load critical N, P losses
        if crop == "winterwheat":
            crop_crit_name = "Wheat"
        elif crop == "maize":
            crop_crit_name = "Maize"
        elif crop == "soybean":
            crop_crit_name = "Soybean"
        elif crop == "mainrice":
            crop_crit_name = "Rice"
        elif crop == "secondrice":
            crop_crit_name = "Rice" 
            
        crit_N_loss_file = os.path.join(crit_loss_dir, f"{crop_crit_name}/{basin}_crit_N_runoff_kg.nc")
        ds_crit_N_loss = xr.open_dataset(crit_N_loss_file)
        Crit_N_Runoff = ds_crit_N_loss["Critical_N_runoff"]

        crit_P_loss_file = os.path.join(crit_loss_dir, f"{crop_crit_name}/{basin}_crit_P_runoff_kg.nc")
        ds_crit_P_loss = xr.open_dataset(crit_P_loss_file)
        Crit_P_Runoff = ds_crit_P_loss["Critical_P_runoff"]

        # Load harvested area
        if crop == "winterwheat":
            cropname = "WHEA"
        elif crop == "maize":
            cropname = "MAIZ"
        elif crop == "soybean":
            cropname = "SOYB"
        elif crop == "mainrice" and basin != "Yangtze":
            cropname = "RICE"
        elif crop == "mainrice" and basin == "Yangtze":
            cropname = "MAINRICE"
        elif crop == "secondrice":
            cropname = "SECONDRICE" 

        Total_HA_file = os.path.join(data_dir, basin, f"Harvest_Area/{cropname}_Harvest_Area_05d_{basin}.nc")
        ds_total_HA = xr.open_dataset(Total_HA_file)
        Total_HA = ds_total_HA["Harvest_Area"]   

        Rainfed_HA_file = os.path.join(data_dir, basin, f"Harvest_Area/{cropname}_Rainfed_Harvest_Area_05d_{basin}.nc")
        ds_Rainfed_HA = xr.open_dataset(Rainfed_HA_file)
        Rainfed_HA = ds_Rainfed_HA["Harvest_Area"]             

        Irrigated_HA_file = os.path.join(data_dir, basin, f"Harvest_Area/{cropname}_Irrigated_Harvest_Area_05d_{basin}.nc")
        ds_Irrigated_HA = xr.open_dataset(Irrigated_HA_file)
        Irrigated_HA = ds_Irrigated_HA["Harvest_Area"]

        # Load irrigation file
        Irrigation_file = os.path.join(data_dir,basin,f"Irrigation_2005-2019/{basin}_{crop}_monthly_Irri_Rate.nc")
        ds_Irrigation = xr.open_dataset(Irrigation_file)
        Irrigation_monthly = ds_Irrigation["Irrigation_Rate"].sel(time=slice("2015-01-01", "2015-12-01"))
        Annual_Irri_rate = Irrigation_monthly.sum(dim= "time")

        sus_Irrigation_file = os.path.join(data_dir,basin,f"Irrigation_2005-2019/{basin}_{crop}_monthly_sus_Irri_Rate.nc")
        ds_sus_Irrigation = xr.open_dataset(sus_Irrigation_file)
        sus_Irrigation_monthly = ds_sus_Irrigation["Irrigation_Rate"].sel(time=slice("2015-01-01", "2015-12-01"))
        Annual_sus_Irri_rate = sus_Irrigation_monthly.sum(dim= "time")

        # Calculate the total losses [kg] and water use
        N_Runoff_total = N_Runoff_Irrigated.fillna(0) * Irrigated_HA.fillna(0) + N_Runoff_Rainfed.fillna(0) * Rainfed_HA.fillna(0)
        Crit_N_total = Crit_N_Runoff.fillna(0)

        P_Runoff_total = P_Runoff_Irrigated.fillna(0) * Irrigated_HA.fillna(0) + P_Runoff_Rainfed.fillna(0) * Rainfed_HA.fillna(0)
        Crit_P_total = Crit_P_Runoff.fillna(0)

        Irri_total = Annual_Irri_rate.fillna(0) * Irrigated_HA.fillna(0) * 10 # Transform to m3 
        sus_Irri_total = Annual_sus_Irri_rate.fillna(0) * Irrigated_HA.fillna(0)  * 10 # Transform to m3 

        # Get the mask where the boundaries have been surpassed
        N_exceedance = mask.copy()
        P_exceedance = mask.copy()
        Irri_exceedance = mask.copy()

        # Irri_exceedance = xr.where(Irri_total > sus_Irri_total, 11, Irri_exceedance) # Commnet this one if it is sustainable irrigated
        N_exceedance = xr.where(N_Runoff_total > Crit_N_total, 11, N_exceedance)
        P_exceedance = xr.where(P_Runoff_total > Crit_P_total, 11, P_exceedance)

        lat = mask["lat"].values
        lon = mask["lon"].values

        ds = xr.Dataset(
            {
               "Basin_mask": (("lat", "lon"), mask.values),
               "N_Runoff": (("lat", "lon"), N_Runoff_total.values),                 # [kg]
               "P_Runoff": (("lat", "lon"), P_Runoff_total.values),                 # [kg]
               "Crit_N_Runoff": (("lat", "lon"), Crit_N_total.values),              # [kg]
               "Crit_P_Runoff": (("lat", "lon"), Crit_P_total.values),              # [kg] 
               "Avg_Yield_Rainfed": (("lat", "lon"), Avg_Yield_Rainfed.values),     # [kg/ha]
               "Avg_Yield_Irrigated": (("lat", "lon"), Avg_Yield_Irrigated.values), # [kg/ha]
               "Total_HA": (("lat", "lon"), Total_HA.values),                       # [ha]
               "Rainfed_HA": (("lat", "lon"), Rainfed_HA.values),                   # [ha]
               "Irrigated_HA": (("lat", "lon"), Irrigated_HA.values),               # [ha]
               "Total_irrigation_amount": (("lat", "lon"), Irri_total.values),      # [m3]
               "Sus_irrigation_amount": (("lat", "lon"), sus_Irri_total.values),    # [m3] 
               "Irri_exceedance": (("lat", "lon"), Irri_exceedance.values),
               "N_exceedance": (("lat", "lon"), N_exceedance.values),
               "P_exceedance": (("lat", "lon"), P_exceedance.values),
            },
            coords = {"lat": lat, "lon": lon}
        )

        ds.attrs["description"] = f"Baseline scenario analysis for {crop} in {basin}"
        output_file = os.path.join(output_dir, f"{basin}_{crop}_summary_baseline.nc")
        ds.to_netcdf(output_file)
        print(f"----> Successfully summarized baseline scenario and saved {output_file}")
