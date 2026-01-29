import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker # New import for axis control

# --- Configuration ---
Studyareas = ["LaPlata", "Indus", "Yangtze", "Rhine"]
InputCrops = ["winterwheat", "maize", "mainrice", "secondrice", "soybean"]
FinalCategories = ["Wheat", "Maize", "Rice", "Soybean"]

input_dir = "/lustre/nobackup/WUR/ESG/zhou111/3_RQ1_Model_Outputs/4_Analysis4Plotting/0_Summary"
output_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/Demo_Plots/MainFigs"

units = ["Irrigation water amount ($m^3$)", "ktons N", "ktons P"]
cat_colors = {"Wheat": "#E3B448", "Rice": "#4793AF", "Maize": "#D85C27",  "Soybean": "#6B8E23"}
crop_mapper = {"winterwheat": "Wheat", "maize": "Maize", "mainrice": "Rice", "secondrice": "Rice", "soybean": "Soybean"}

# Adjust figsize: width reduced to 12 for skinnier subplots
fig, axes = plt.subplots(4, 3, figsize=(12, 22))
plt.subplots_adjust(hspace=0.35, wspace=0.45, top=0.92)

for b_idx, basin in enumerate(Studyareas):
    basin_data = {m: {cat: 0.0 for cat in FinalCategories} for m in 
                  ["Irri", "Sus_Irri", "N_Runoff", "Crit_N", "P_Runoff", "Crit_P"]}

    for crop in InputCrops:
        file_path = os.path.join(input_dir, f"{basin}_{crop}_summary_baseline.nc")
        if not os.path.exists(file_path): continue
        target_cat = crop_mapper[crop]
        with xr.open_dataset(file_path) as ds:
            mask, i_m, t_m = ds['Basin_mask'], ds['Irrigated_HA'] > 2500, ds['Total_HA'] > 2500
            basin_data["Irri"][target_cat] += (ds['Total_irrigation_amount'].where(i_m) * mask).sum().item()
            basin_data["Sus_Irri"][target_cat] += (ds['Sus_irrigation_amount'].where(i_m) * mask).sum().item()
            basin_data["N_Runoff"][target_cat] += (0.000001 * ds['N_Runoff'].where(t_m) * mask).sum().item()
            basin_data["Crit_N"][target_cat] += (0.000001 *ds['Crit_N_Runoff'].where(t_m) * mask).sum().item()
            basin_data["P_Runoff"][target_cat] += (0.000001 *ds['P_Runoff'].where(t_m) * mask).sum().item()
            basin_data["Crit_P"][target_cat] += (0.000001 *ds['Crit_P_Runoff'].where(t_m) * mask).sum().item()

    metric_cols = [(basin_data["Irri"], basin_data["Sus_Irri"]),
                   (basin_data["N_Runoff"], basin_data["Crit_N"]),
                   (basin_data["P_Runoff"], basin_data["Crit_P"])]


    for m_idx, (base_dict, crit_dict) in enumerate(metric_cols):
        ax = axes[b_idx, m_idx]
        b_base, b_crit = 0.0, 0.0
        bar_width = 0.4 
        
        for cat in FinalCategories:
            v_base = base_dict[cat]
            v_crit = crit_dict[cat]
            ax.bar(0, v_base, bottom=b_base, width=bar_width, color=cat_colors[cat], edgecolor='white', linewidth=0.5)
            ax.bar(1, v_crit, bottom=b_crit, width=bar_width, color=cat_colors[cat], alpha=1.0, edgecolor='white', linewidth=0.5)
            b_base += v_base
            b_crit += v_crit

        # Formatting Cleanup
        ax.set_title("")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel(f"{units[m_idx]}", fontsize=15)
        
        # --- "Nice Number" Limit Logic ---
        max_h = max(b_base, b_crit)
        if max_h > 0:
            # Determine the order of magnitude
            mag = 10**np.floor(np.log10(max_h))
            for nice_step in [1, 2, 5, 10]:
                if max_h <= nice_step * mag:
                    ymax = nice_step * mag
                    break
            else:
                ymax = np.ceil(max_h/mag) * mag
                
            ax.set_ylim(0, ymax)

            # --- Forced Integer Ticks ---
            if m_idx == 0:  # Water Use (Scientific)
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
                # MaxNLocator(integer=True) prevents 0.5, 1.5, etc.
                ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, integer=True))
            else:           # N & P (Normal Values)
                ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.1f}'))
                ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, integer=True))
            
            ax.tick_params(axis='y', labelsize=13)


        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Baseline", "Critical"], fontsize=15)
        ax.set_xlim(-0.6, 1.6)

    # --- Print Detailed Breakdown for each Crop ---
    print(f"\n{'='*75}")
    print(f" BASIN: {basin.upper()}")
    print(f"{'='*75}")
    print(f"{'Crop':<12} | {'Water Use (m3)':>18} | {'N Runoff (ktons)':>15} | {'P Runoff (ktons)':>15}")
    print(f"{'-'*75}")
    
    for cat in FinalCategories:
        w = basin_data["Irri"][cat]
        n = basin_data["N_Runoff"][cat]
        p = basin_data["P_Runoff"][cat]
        # :,.0f adds commas and removes decimals for a "normal" look
        print(f"{cat:<12} | {w:18,.2f} | {n:15,.2f} | {p:15,.2f}")
    print(f"{'-'*75}\n")

# Global legend
legend_elements = [Line2D([0], [0], color=cat_colors[cat], lw=8, label=cat) for cat in FinalCategories]
# Create the legend at the bottom
fig.legend(
    handles=legend_elements, 
    loc='lower center',       # Align to the lower center
    bbox_to_anchor=(0.5, 0.05), # x=0.5 (middle), y=0.05 (just above the bottom edge)
    ncol=4,                   # Keep them in one row
    frameon=False, 
    fontsize=16
)

plt.savefig(os.path.join(output_dir, "Basin_Slim_Integral_Final.png"), bbox_inches='tight', dpi=300)

