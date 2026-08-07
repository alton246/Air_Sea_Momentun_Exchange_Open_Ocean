import pandas as pd
import numpy as np
import xarray as xr
import cmaps
import matplotlib.pyplot as plt

from helpers import *

PATH = '/home/disk/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/notebooks/' # Path to data
PNG = '/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/' # Path to save figures
file = 'storm_relative_qualit_controlled_data.*.nc' # File pattern

data = xr.open_mfdataset(PATH + file, combine='by_coords') # Open multiple netCDF files as a single xarray dataset

r = np.sqrt(data.x**2 + data.y**2) # Compute radial distance
dist = r.values.flatten() # Flatten distance array

awo_stress_ratio_bins = np.linspace(0, 1.6, 50) # Define stress ratio bins
dist_bins = np.linspace(20, 200, 100) # Define distance bins

awo_stress_ratio_mean_comb = [] # Initialize list to store mean CTL stress ratios
awo_ws_stress_ratio_mean_comb = [] # Initialize list to store EXP stress ratios


mid_periods = pd.date_range(start='2010-08-30T15:00:00', end='2010-09-01T09:00:00', 
                            freq='6h') # Define mid-periods for analysis
n = len(mid_periods) # Number of mid-periods

#Figure Settings
gridsize = (4, 1)

colors=['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet', 'k']
xlim = (0, 200)
ylim = (0.75, 1.10)
fontsize=4.5
labelsize=4.5
labelpad=0.25
linewidth=0.5
alpha=0.5
s=0.1
lw=0.8
zorder=50
xpos, ypos=0.044, 0.84

awo_labels = ['(a)', '(c)', '(e)',  '(g)']
awo_ws_labels = ['(b)', '(d)', '(f)', '(h)']
fig = plt.figure(figsize=(2, 4))

print(data['time'].values)

for i in range(len(data['time'][1:][::2])):

    print(data['time'][1:][::2][i].values)

    # Get the stress fields for the current time step and flatten them
    awo_ocn_stress_ctl = data['awo_ocn_stress'][1:][::2][i].values.flatten()
    awo_ocn_stress_exp = data['awo_ws_ocn_stress'][1:][::2][i].values.flatten()

    awo_atm_stress_ctl = data['awo_atm_stress'][1:][::2][i].values.flatten()
    awo_atm_stress_exp = data['awo_ws_atm_stress'][1:][::2][i].values.flatten()

    #Locate Radius of Maximum Wind
    awo_max_dist = get_max_wind_distance(dist, data['wspd_awo_qc'][1:][::2][i].values.flatten())
    awo_ws_max_dist = get_max_wind_distance(dist, data['wspd_awo_ws_qc'][1:][::2][i].values.flatten())

    awo_dist_bins = np.linspace(awo_max_dist, 200, 100)
    awo_ws_dist_bins = np.linspace(awo_ws_max_dist, 200, 100)

    
    awo_stress_ratio = awo_ocn_stress_ctl / awo_atm_stress_ctl #CTL
    awo_ws_stress_ratio = awo_ocn_stress_exp / awo_atm_stress_exp #EXP

    # Bin the distances and compute average ratio at each radius
    awo_stress_ratio_means, awo_stress_ratio_std, awo_hist, awo_xedges, awo_yedges = compute_stress_ratio_means(dist, awo_stress_ratio, dist_bins=awo_dist_bins, ratio_bins=awo_stress_ratio_bins)
    awo_ws_stress_ratio_means, awo_ws_stress_ratio_std, awo_ws_hist, awo_ws_xedges, awo_ws_yedges = compute_stress_ratio_means(dist, awo_ws_stress_ratio, dist_bins=awo_ws_dist_bins, ratio_bins=awo_stress_ratio_bins)

    awo_std_plus_mean = np.array(awo_stress_ratio_std) + np.array(awo_stress_ratio_means)
    awo_std_minus_mean = np.array(awo_stress_ratio_means) - np.array(awo_stress_ratio_std)
    awo_ws_std_plus_mean = np.array(awo_ws_stress_ratio_std) + np.array(awo_ws_stress_ratio_means)
    awo_ws_std_minus_mean = np.array(awo_ws_stress_ratio_means) - np.array(awo_ws_stress_ratio_std)

    # Create subplots for each time step: CTL (left), EXP (right)
    row = i % gridsize[0]  # Ensure row index is within grid size
    ax1 = plt.subplot2grid(gridsize, (row, 0), colspan=1, rowspan=1)
    # ax2 = plt.subplot2grid(gridsize, (row, 1), colspan=1, rowspan=1)

    # CTL (AWO) on the left
    ax1.scatter(dist, awo_stress_ratio, s=s, color='k', alpha=0.5)
    ax1.scatter(dist, awo_ws_stress_ratio, s=s, color='blue', alpha=0.5)
    ax1.plot(awo_ws_xedges, awo_ws_stress_ratio_means, color='blue', label='EXP mean', zorder=zorder, lw=lw)
    ax1.plot(awo_xedges, awo_stress_ratio_means, color='orange', lw=lw, label='CTL mean')
    ax1.plot(awo_xedges, awo_std_plus_mean, color='orange', lw=0.5, label=r'$\pm$ CTL std', linestyle='--', zorder=zorder)
    ax1.plot(awo_xedges, awo_std_minus_mean, color='orange', lw=0.5, linestyle='--', zorder=zorder)
    ax1.plot(awo_ws_xedges, awo_ws_std_plus_mean, color='blue', lw=0.5, label=r'$\pm$ EXP std', linestyle='--', zorder=zorder)
    ax1.plot(awo_ws_xedges, awo_ws_std_minus_mean, color='blue', lw=0.5, linestyle='--', zorder=zorder)
    ax1.hlines(y=1, xmin=xlim[0], xmax=xlim[1], colors='red', linestyles='--', 
               label='$\\tau_{ocn} = \\tau_{atm}$', lw=lw, zorder=200)
    ax1.vlines(x=awo_max_dist, ymin=0, ymax=ylim[1], colors='k', linestyles='-', lw=lw, label='CTL RMW')
    ax1.vlines(x=awo_ws_max_dist, ymin=0, ymax=ylim[1], colors='green', linestyles='-', lw=lw, label='EXP RMW')

    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_xlabel('Distance from Storm Center (km)', fontsize=fontsize, labelpad=labelpad)
    ax1.set_ylabel('Stress Ratio ($\\tau_{ocn} / \\tau_{atm}$)', fontsize=fontsize, labelpad=labelpad)
    ax1.set_title(mid_periods[1:][::2][i].strftime('%Y-%m-%d %H'), fontsize=fontsize, pad=labelpad)
    ax1.tick_params(axis='both', which='major', labelsize=labelsize)
    ax1.grid(linestyle='--', linewidth=linewidth, alpha=alpha)
    # ax1.legend(fontsize=3, loc='best', framealpha=0.8, shadow=True)
    add_corner_label(ax1, xpos, ypos, awo_labels[i], fontsize=4)


    handles1, labels1 = ax1.get_legend_handles_labels()
    # handles2, labels2 = ax2.get_legend_handles_labels()
    fig.legend(handles1, labels1, loc='upper center', ncol=3, fontsize=4, frameon=True, 
               bbox_to_anchor=(0.5, 1.1), shadow=True, framealpha=0.8)

fig.tight_layout(pad=0.01, w_pad=1, h_pad=0.5)

plt.savefig(PNG + 'fig6_mom_trans_fact_scatter_ESS.png', dpi=500, bbox_inches='tight',
                facecolor='w', transparent=False)