# This script is used to plot hte P pool [mg/kg] v.s. P runoff [kg/ha] for 4 basins

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

# Configuration
# Studyareas = ["LaPlata", "Indus"] # ["LaPlata", "Indus", "Yangtze", "Rhine"]
# ymin, ymax = 20, 120 # Minimum and Maximum P pool size

# Studyareas = ["Rhine"] # ["LaPlata", "Indus", "Yangtze", "Rhine"]
# ymin, ymax = 500, 650 # Minimum and Maximum P pool size

# Studyareas = ["Yangtze"] 
# ymin, ymax = 150, 350 # Minimum and Maximum P pool size

Studyareas = ["Indus"] 
ymin, ymax = 0, 100 # Minimum and Maximum P pool size

CropTypes = ["winterwheat"] # ["winterwheat", "mainrice", "secondrice", "maize", "soybean"]

file_dir = "/lustre/nobackup/WUR/ESG/zhou111/4_RQ1_Analysis_Results/V2_Demo_Plots/Fig5/5c"


def linear_zero_intercept(x, m):
    return m * x



for basin in Studyareas:
    for crop in CropTypes:
        file_path = os.path.join(file_dir, f"{basin}_{crop}_P_analysis.csv")

        if not os.path.exists(file_path):
            continue

        df = pd.read_csv(file_path)
        df['DisplayYear'] = df['Year'] + 1
        target_years = [2010, 2015, 2020]

        scenarios = {
            "Baseline": {"pool": "Baseline_P_Pool", "runoff": "Baseline_P_Runoff", "color":"orange", "label": "Sustainable irrigation scenario"},
            "Red_prop": {"pool": "Red_prop_P_Pool", "runoff": "Red_prop_P_Runoff", "color": "teal", "label": "Fertilizer reduction scenario"},
        }

        # 1) Use gridspec_kw to make Subplot 1 (Time Series) twice as wide as Subplot 2 (Scatter)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), 
                                       gridspec_kw={'width_ratios': [2, 1]})

        # --- SUBPLOT 1: Time Series Year vs P Pool (Wider Plot) ---
        for key, config in scenarios.items():
            ax1.plot(df['DisplayYear'], df[config["pool"]], color=config["color"], 
                     linewidth=4, label=config["label"], marker='o', markersize=12, 
                     markeredgecolor='white', markeredgewidth=1.5)
            
        ax1.set_ylim(ymin, ymax)
        ax1.set_xlabel("Year", fontsize=18)
        ax1.set_ylabel("Basin Average P Pool [mg P/kg]", fontsize=18)
        ax1.set_xticks(df['DisplayYear'])
        ax1.tick_params(axis='both', labelsize=16)
        ax1.grid(True, linestyle='--', alpha=0.5)
        ax1.legend(fontsize=18, loc='upper left')

        # --- SUBPLOT 2: Scatter P Runoff (X) vs P Pool (Y) (Swapped Axes) ---
        all_pool = []
        all_runoff = []

        for key, config in scenarios.items():
            # X is now Runoff, Y is now Pool
            x_data = df[config["runoff"]]
            y_data = df[config["pool"]]
            
            all_runoff.extend(x_data)
            all_pool.extend(y_data)
            
            ax2.scatter(x_data, y_data, color=config["color"], alpha=0.8, edgecolors='w', s=200, zorder=3)
            
            # Labels for specific years
            for i, row in df.iterrows():
                if row['DisplayYear'] in target_years:
                    ax2.text(row[config["runoff"]], row[config["pool"]], f" {int(row['DisplayYear'])}", 
                             fontsize=16, color="black")

        # Combined Trendline for Subplot 2 (forcing intercept 0)
        # Note: We fit y (Pool) as a function of x (Runoff)
        popt, _ = curve_fit(linear_zero_intercept, all_runoff, all_pool)
        m_slope = popt[0]
        
        x_fit_range = np.linspace(min(all_runoff)*0.8, max(all_runoff)*1.2, 100)
        ax2.set_ylim(ymin, ymax)
        ax2.plot(x_fit_range, linear_zero_intercept(x_fit_range, m_slope), 
                 color='black', linestyle=':', linewidth=2, alpha=0.7, 
                 label=f'Slope: {m_slope:.2f}')
        ax2.set_xlim(0, 1.5)

        ax2.set_xlabel("P Runoff [kg/ha]", fontsize=18)
        ax2.set_ylabel("P Pool [mg P/kg]", fontsize=18)
        ax2.tick_params(axis='both', labelsize=14)
        ax2.grid(True, linestyle='--', alpha=0.5)

        # Final adjustments
        plt.tight_layout()
        save_name = os.path.join(file_dir, f"{basin}_{crop}_final_v3.png")
        plt.savefig(save_name, dpi=300)
        plt.show()

        print(f"Updated analysis plot saved for {basin} {crop}.")