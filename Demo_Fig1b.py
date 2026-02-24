import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator

# --- Settings & Directories ---
Studyarea = ["LaPlata", "Rhine", "Indus", "Yangtze"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

model_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

fig_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig1/1b"
os.makedirs(fig_output_dir, exist_ok=True)

Variables = {
    "Water": {"label": "Irrigation Exceedance (m3)", "cmap": "Blues"},
    "N": {"label": "N Runoff Exceedance (ktons)", "cmap": "RdPu"},
    "P": {"label": "P Runoff Exceedance (ktons)", "cmap": "YlOrRd"}
}

def GetCropWNP(file_name):
    ds = xr.open_dataset(file_name)
    mask = ds["Basin_mask"].where(ds["Total_HA"] > 2500, 0)
    
    Total_HA = ds["Total_HA"] * mask
    Total_Irri = ds["Total_irrigation_amount"] * mask
    Sus_Irri   = ds["Sus_irrigation_amount"] * mask
    N_runoff   = 0.000001 * ds["N_Runoff"] * mask
    P_runoff   = 0.000001 * ds["P_Runoff"] * mask
    Crit_N     = 0.000001 * ds["Crit_N_Runoff"] * mask
    Crit_P     = 0.000001 * ds["Crit_P_Runoff"] * mask
    
    return Total_HA, Total_Irri, Sus_Irri, N_runoff, P_runoff, Crit_N, Crit_P

def plot_exceedance_map(var_key, data_dict, basins, output_path):
    # 1. Calculate spans for normalization (consistent height logic)
    max_lat_span = 0
    basin_info = {}
    for basin in basins:
        gdf = gpd.read_file(os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp"))
        minx, miny, maxx, maxy = gdf.total_bounds
        basin_info[basin] = {'bounds': (minx, miny, maxx, maxy), 
                             'lat_span': maxy - miny, 'lon_span': maxx - minx}
        if (maxy - miny) > max_lat_span: max_lat_span = maxy - miny

    # 2. Layout Setup
    ratios = [basin_info[b]['lon_span'] for b in basins]
    fig = plt.figure(figsize=(26, 8)) # Slightly wider to accommodate right colorbar
    gs = gridspec.GridSpec(1, 4, width_ratios=ratios, figure=fig, wspace=0.5)
    
    # 3. Colormap Logic (Categorical Green for Sustainable)
    base_cmap = plt.get_cmap(Variables[var_key]["cmap"]).copy()
    base_cmap.set_under("#dffad9ea") # Your light green categorical color
    
    all_values = [da.values.flatten() for da in data_dict.values() if da is not None]
    flat_values = np.concatenate(all_values)
    vmax = np.nanmax(flat_values)
    
    # 1e-10 ensures 0 and negatives fall into 'set_under'
    norm = mcolors.Normalize(vmin=1e-10, vmax=vmax)

    for i, basin in enumerate(basins):
        ax = fig.add_subplot(gs[i], projection=ccrs.PlateCarree())
        ax.spines['geo'].set_visible(False)
        
        # 4. Background Features
        ax.add_feature(cfeature.LAND, facecolor="#f5f5f5f2")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor="gray")
        
        # 5. Low Runoff Mask
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
        if os.path.exists(low_runoff_path):
            with xr.open_dataset(low_runoff_path) as ds_lr:
                lr_mask = ds_lr["Low_Runoff"] 
            lr_mask.plot.imshow(ax=ax, transform=ccrs.PlateCarree(), 
                                cmap=mcolors.ListedColormap(["#C4C4C4FF"]), 
                                add_colorbar=False, zorder=2)

        # 6. Plot Main Data
        data = data_dict[basin]
        im = data.plot.imshow(
            ax=ax, transform=ccrs.PlateCarree(),
            norm=norm, cmap=base_cmap,
            add_colorbar=False, zorder=4
        )
        
        # 7. Basin Boundary
        gdf = gpd.read_file(os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp"))
        gdf.boundary.plot(ax=ax, color='black', linewidth=1.5, zorder=5)

        # 8. Set View Extent
        minx, miny, maxx, maxy = basin_info[basin]['bounds']
        center_y = (maxy + miny) / 2
        view_miny, view_maxy = center_y - (max_lat_span/2) - 2, center_y + (max_lat_span/2) + 2
        ax.set_extent([minx - 2, maxx + 2, view_miny, view_maxy], crs=ccrs.PlateCarree())
        
        # 9. Ticks/Coordinates
        gl = ax.gridlines(draw_labels=True, linewidth=0, color='none')
        gl.top_labels = False; gl.right_labels = False
        gl.xlocator = MultipleLocator(5); gl.ylocator = MultipleLocator(5)
        gl.xlabel_style = {'size': 15, 'color': 'black'}
        gl.ylabel_style = {'size': 15, 'color': 'black'}

    # 10. Vertical Colorbar on the Right
    # [left, bottom, width, height]
    # left=0.92 puts it outside the subplots area
    cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5]) 
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar.outline.set_visible(False)
    
    cbar.set_label(Variables[var_key]["label"], fontsize=20, labelpad=20)
    cbar.ax.tick_params(labelsize=15)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()

# --- Main Processing Loop ---
results = { "Water": {}, "N": {}, "P": {} }

for basin in Studyarea:
    print(f"Processing Basin: {basin}")

    low_runoff_path = os.path.join(data_dir, "2_StudyArea", basin, "low_runoff_mask.nc")
    with xr.open_dataset(low_runoff_path) as ds_lr:
        low_runoff = ds_lr["Low_Runoff"]
    mask_not_low_runoff = xr.where(low_runoff.isnull(), 1, np.nan)

    basin_W = 0; basin_W_crit = 0
    basin_N = 0; basin_N_crit = 0
    basin_P = 0; basin_P_crit = 0
    basin_total_HA = 0
    
    found_data = False
    for crop in Croptypes:
        file_name = os.path.join(model_output_dir, f"{basin}_{crop}_summary.nc")
        if not os.path.exists(file_name):
            continue 
        
        HA, W, Wc, N, P, Nc, Pc = GetCropWNP(file_name)
        
        basin_W += W.fillna(0)
        basin_W_crit += Wc.fillna(0)
        basin_N += N.fillna(0)
        basin_N_crit += Nc.fillna(0)
        basin_P += P.fillna(0)
        basin_P_crit += Pc.fillna(0)
        basin_total_HA += HA.fillna(0)
        found_data = True
    
    mask_low_HA = xr.where(basin_total_HA > 2500, 1, np.nan)
    results["Water"][basin] = (basin_W - basin_W_crit) * mask_not_low_runoff * mask_low_HA
    results["N"][basin] = (basin_N - basin_N_crit) * mask_not_low_runoff * mask_low_HA
    results["P"][basin] = (basin_P - basin_P_crit) * mask_not_low_runoff * mask_low_HA

for var_key in Variables:
    out_file = os.path.join(fig_output_dir, f"Exceedance_{var_key}.png")
    plot_exceedance_map(var_key, results[var_key], Studyarea, out_file)
    print(f"Saved figure for {var_key}")