import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

Studyarea =  ["Indus"] # ["Indus", "LaPlata", "Yangtze", "Rhine"]
Croptypes =  ["winterwheat"] # ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data/2_StudyArea"
input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/1_Summary"

for basin in Studyarea:
    basin_shp_file = os.path.join(data_dir, f"2_shp_StudayArea/{basin}/{basin}.shp")
    river_shp_file = os.path.join(data_dir, f"2_shp_StudayArea/{basin}/{basin}_River.shp")
    for crop in Croptypes:
        if crop == "secondrice": 
            input_file1 = os.path.join(input_dir, f"{basin}_mainrice_baseline_summary")
            input_file2 = os.path.join(input_dir, f"{basin}_secondrice_baseline_summary")
       
        else:
            input_file = os.path.join(input_dir, f"{basin}_{crop}_baseline_summary")