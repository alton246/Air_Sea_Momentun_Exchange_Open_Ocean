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

for i in range(len(data['time'])):

    # Get the stress fields for the current time step and flatten them
    awo_ocn_stress_ctl = data['awo_ocn_stress'][i].values.flatten()
    awo_ocn_stress_exp = data['awo_ws_ocn_stress'][i].values.flatten()

    awo_atm_stress_ctl = data['awo_atm_stress'][i].values.flatten()
    awo_atm_stress_exp = data['awo_ws_atm_stress'][i].values.flatten()

    # r = np.sqrt(data.x**2 + data.y**2)
    # Compute the stress ratios
    awo_stress_ratio = awo_ocn_stress_ctl / awo_atm_stress_ctl #CTL
    awo_ws_stress_ratio = awo_ocn_stress_exp / awo_atm_stress_exp #EXP

    # Bin the distances and compute average ratio at each radius
    awo_stress_ratio_means, awo_stress_ratio_std, awo_hist, awo_xedges, awo_yedges = compute_stress_ratio_means(dist, awo_stress_ratio, dist_bins=dist_bins, ratio_bins=awo_stress_ratio_bins)
    awo_ws_stress_ratio_means, awo_ws_stress_ratio_std, awo_ws_hist, awo_ws_xedges, awo_ws_yedges = compute_stress_ratio_means(dist, awo_ws_stress_ratio, dist_bins=dist_bins, ratio_bins=awo_stress_ratio_bins)

    awo_stress_ratio_mean_comb.append([float(x) for x in awo_stress_ratio_means]) # Append CTL means to list
    awo_ws_stress_ratio_mean_comb.append([float(x) for x in awo_ws_stress_ratio_means]) # Append EXP means to list

    awo_xedges_list = [awo_xedges.copy() for _ in range(n)]
    awo_ws_xedges_list = [awo_ws_xedges.copy() for _ in range(n)]


#Figure Settings
gridsize = (1, 2)

colors=['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet', 'k']
xlim = (0, 200)
ylim = (0.80, 1.05)
fontsize=8
labelsize=6
labelpad=0.25
linewidth=0.5
alpha=0.5
lw=0.8
zorder=50
xpos, ypos=0.082, 0.92

#Plotting Figure
fig = plt.figure(figsize=(4, 2))
#AWO
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)

for idx in range(len(awo_xedges_list)):
    ax1.plot(awo_xedges_list[idx], awo_stress_ratio_mean_comb[idx], color=colors[idx], lw=lw, 
             label=mid_periods[idx].strftime('%Y-%m-%d %H:%M')) # Plot CTL stress ratio means for eaxh period
ax1.hlines(y=1, xmin=xlim[0], xmax=xlim[1], colors='red', linestyles='--', 
           label='$\\tau_{ocn} = \\tau_{atm}$', lw=lw) # Reference line at y=1

ax1.set_xlim(xlim) # Set x-axis limits
ax1.set_ylim(ylim) # Set y-axis limits

ax1.set_xlabel('Distance from Storm Center (km)', fontsize=fontsize, labelpad=labelpad) # Set x-axis label
ax1.set_ylabel('Stress Ratio ($\\tau_{ocn} / \\tau_{atm}$)', fontsize=fontsize, labelpad=labelpad) # Set y-axis label
ax1.set_title('CTL', fontsize=fontsize, pad=labelpad) # Set title
ax1.tick_params(axis='both', which='major', labelsize=labelsize) # Set tick parameters
ax1.grid(linestyle='--', linewidth=linewidth, alpha=alpha) # Add grid
add_corner_label(ax1, xpos, ypos, '(a)', fontsize=4)

ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)

for idx in range(len(awo_xedges_list)):
    ax2.plot(awo_ws_xedges_list[idx], awo_ws_stress_ratio_mean_comb[idx], color=colors[idx], 
             lw=lw, zorder=zorder, label=mid_periods[idx].strftime('%Y-%m-%d %H:%M')) # Plot EXP stress ratio means for each period

ax2.hlines(y=1, xmin=xlim[0], xmax=xlim[1], colors='red', linestyles='--', 
           label='$\\tau_{ocn} = \\tau_{atm}$', lw=lw) # Reference line at y=1

ax2.set_xlim(xlim) # Set x-axis limits
ax2.set_ylim(ylim) # Set y-axis limits

ax2.set_xlabel('Distance from Storm Center (km)', fontsize=fontsize, labelpad=labelpad) # Set x-axis label
ax2.set_ylabel('Stress Ratio ($\\tau_{ocn} / \\tau_{atm}$)', fontsize=fontsize, labelpad=labelpad) # Set y-axis label
ax2.set_title('EXP', fontsize=fontsize, pad=labelpad) # Set title
ax2.tick_params(axis='both', which='major', labelsize=labelsize) # Set tick parameters
ax2.grid(linestyle='--', linewidth=linewidth, alpha=alpha) # Add grid
ax2.legend(fontsize=3.5, loc='lower right', frameon=True, shadow=True) # Add legend
add_corner_label(ax2, xpos, ypos, '(b)', fontsize=4)
# ax2.set_aspect('equal') # Set aspect ratio to equal

fig.tight_layout(pad=0.01, w_pad=0.5, h_pad=-4)


plt.savefig(PNG + 'fig_ocn_stress_evolution_JAMES.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)