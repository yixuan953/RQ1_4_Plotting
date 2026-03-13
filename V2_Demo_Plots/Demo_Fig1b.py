import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import cartopy.geodesic as cgeo
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LinearSegmentedColormap

# --- Settings & Directories ---
Studyarea = ["LaPlata", "Rhine", "Indus", "Yangtze"]
Croptypes = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]

model_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary/1_Baseline"
data_dir = "/lustre/nobackup/WUR/ESG/zhou111/2_RQ1_Data"

fig_output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig1/1b"
os.makedirs(fig_output_dir, exist_ok=True)

p_colors = ["#fff6e9f1", "#feb728"]
custom_phosphorus_cmap = LinearSegmentedColormap.from_list("custom_p", p_colors, N=256)

Variables = {
    "Water": {"label": "Irrigation Exceedance (m3)", "cmap": "Blues"},
    "N": {"label": "N Runoff Exceedance (ktons)", "cmap": "RdPu"},
    "P": {"label": "P Runoff Exceedance (ktons)", "cmap": custom_phosphorus_cmap}
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

def add_scale_bar(ax, basin, length_km=None, location_y=-0.08):
    """
    Places a scale bar dynamically based on basin size.
    """
    lon0, lon1, lat0, lat1 = ax.get_extent()
    
    # Determine appropriate scale bar length if not provided
    if length_km is None:
        span_lon = lon1 - lon0
        length_km = 500 if span_lon > 10 else 100
        if span_lon < 2: length_km = 20 # Small basins like Rhine
        
    center_lat = (lat0 + lat1) / 2
    geod = cgeo.Geodesic()
    # Calculate how many meters are in 1 degree of longitude at this latitude
    dist_1deg = geod.inverse((lon0, center_lat), (lon0 + 1, center_lat))[0, 0]
    bar_width_deg = (length_km * 1000) / dist_1deg
    
    # Logic for side placement (matching your reference)
    if basin in ["Yangtze"]:
        x_start = lon0 + (lon1 - lon0) * 0.1
        x_end = x_start + bar_width_deg
        ha_text = 'left'
        text_x = x_start
    else: 
        x_end = lon1 - (lon1 - lon0) * 0.1
        x_start = x_end - bar_width_deg
        ha_text = 'right'
        text_x = x_end

    y_pos = lat0 + (lat1 - lat0) * location_y
    
    # Draw scale bar
    ax.plot([x_start, x_end], [y_pos, y_pos], transform=ccrs.PlateCarree(),
            color='black', linewidth=3, zorder=10, clip_on=False)
    
    # Add text
    ax.text(text_x, y_pos - (lat1 - lat0) * 0.02, f'{length_km} km', 
            transform=ccrs.PlateCarree(), ha=ha_text, va='top', 
            fontsize=25, clip_on=False)

def plot_exceedance_map(var_key, data_dict, basins, output_path):
    # 1. Calculate the aspect ratio (width/height) for every basin
    basin_info = []
    total_relative_width = 0
    padding = 0.05  # Space between subplots
    
    for basin in basins:
        gdf = gpd.read_file(os.path.join(data_dir, f"2_shp_StudyArea/{basin}/{basin}.shp"))
    
        minx, miny, maxx, maxy = gdf.total_bounds
        if basin == "Rhine":
            maxy += 0.75 
            miny -= 0.75 
        
        # aspect = width / height
        aspect = (maxx - minx) / (maxy - miny)
        
        basin_info.append({
            'name': basin,
            'bounds': (minx, maxx, miny, maxy),
            'aspect': aspect,
            'gdf': gdf
        })
        total_relative_width += aspect

    # 2. Create Figure
    # We scale the figure width based on the sum of the basin aspects
    fig_width = 26
    fig = plt.figure(figsize=(fig_width, 8))
    
    # 3. Manually place axes to ensure identical height
    # Available width for all basins (1.0 minus total padding and colorbar space)
    available_width = 0.85 - (padding * (len(basins) - 1))
    current_x = 0.05 # Left margin
    
    ims = []
    for info in basin_info:
        # Calculate the width this basin needs to maintain its shape at fixed height
        rel_width = (info['aspect'] / total_relative_width) * available_width
        
        # [left, bottom, width, height] -> height is ALWAYS 0.7
        ax_rect = [current_x, 0.15, rel_width, 0.7]
        ax = fig.add_axes(ax_rect, projection=ccrs.PlateCarree())
        
        current_x += rel_width + padding # Move to the next slot
        
        # 4. Formatting
        ax.axis('off')
        ax.spines['geo'].set_visible(False)
        minx, maxx, miny, maxy = info['bounds']
        ax.set_extent([minx, maxx, miny, maxy], crs=ccrs.PlateCarree())

        # 5. Plotting
        # Low Runoff Mask
        low_runoff_path = os.path.join(data_dir, "2_StudyArea", info['name'], "low_runoff_mask.nc")
        if os.path.exists(low_runoff_path):
            with xr.open_dataset(low_runoff_path) as ds_lr:
                ds_lr["Low_Runoff"].plot.imshow(ax=ax, transform=ccrs.PlateCarree(), 
                                               cmap=mcolors.ListedColormap(["#C4C4C4"]), 
                                               add_colorbar=False, zorder=2)

        # Main Data
        vmax = np.nanmax([da.max() for da in data_dict.values() if da is not None]) if var_key != "Water" else 5e8
        norm = mcolors.Normalize(vmin=1e-10, vmax=vmax)
        base_cmap = plt.get_cmap(Variables[var_key]["cmap"]).copy()
        base_cmap.set_under("#dffad9ea")

        im = data_dict[info['name']].plot.imshow(ax=ax, transform=ccrs.PlateCarree(), 
                                                norm=norm, cmap=base_cmap, 
                                                add_colorbar=False, zorder=4)
        ims.append(im)
        info['gdf'].boundary.plot(ax=ax, color='black', linewidth=1.5, zorder=5)

        # Scale Bar
        add_scale_bar(ax, info['name'])

    # 6. Colorbar (Fixed position on the far right)
    cbar_ax = fig.add_axes([0.92, 0.25, 0.012, 0.5]) 
    cbar = fig.colorbar(ims[-1], cax=cbar_ax, orientation='vertical', extend='both')
    cbar.outline.set_visible(False)
    cbar.set_label(Variables[var_key]["label"], fontsize=20, labelpad=20)

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