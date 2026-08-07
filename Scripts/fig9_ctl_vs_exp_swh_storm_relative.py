import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime 
# import mompy as mom
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import matplotlib
import xarray as xr
import cmaps
import pylab
import os
import glob
import numpy.ma as ma

from datetime import datetime,timedelta
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rcParams
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy import interpolate
# from mompy import *
from helpers import *
from matplotlib.patches import Wedge

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt


DATA_PATH = '/home/disk/orca/adaley17/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/notebooks/'
SUM_PATH = '/home/disk/orca/adaley17/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/'
PNG = '/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/' # Path to save figures
data = 'storm_relative_qualit_controlled_data.*.nc'

radius_df = pd.read_csv(DATA_PATH + 'radius_df.csv') 

print(radius_df.head())

qual_data = xr.open_mfdataset(DATA_PATH + data, combine='by_coords')
ctl_mod_radius_data = extract_storm_radius_data(SUM_PATH + 'earl_5.log_awo.csv')
exp_mod_radius_data = extract_storm_radius_data(SUM_PATH + 'earl_5.log_awo-ws.csv')

os.chdir(DATA_PATH)

file_names = sorted(glob.glob(data))


mid_periods = pd.date_range(start='2010-08-30T15:00:00', end='2010-09-01T09:00:00', 
                            freq='6h') # Define mid-periods for analysis


xmin = -250
xmax = 250
ymin = xmin
ymax = xmax

ticks=np.arange(-250, 300, 50)

xloc, yloc = 0.05, 0.88
#Figure Settings
fontsize=10
labelsize=10
gridsize =(4,3)
s=10
lw=1
y_pos = 0.87
x_pos = 0.090
# zorder=200
# size=50

#Colorbar Settings
wspd_levels = np.arange(5,57.5,2.5)
CMAPS = cmaps.MPL_gist_rainbow_r
shrink=0.45
aspect=18
width=0.5
length=4

labelpad = 0.25

#Color Maps
wspd_cmaps = cmaps.cmp_haxby_r
swh_cmaps = cmaps.NMCRef
diff_cmap = cmaps.MPL_bwr
sst_cmap = cmaps.MPL_jet
sfc_cmap=cmaps.BkBlAqGrYeOrReViWh200_r
stress_cmap = cmaps.precip3_16lev

#Levels
wspd_diff_levels = np.arange(-4,5,1)
swh_diff_levels = np.arange(-0.7,0.8,0.1) 
sst_diff_levels = np.arange(-0.40,0.45, 0.05)
sfc_diff_levels = np.arange(-0.5,0.55,0.05)
wspd_levels = np.arange(5,47.5,2.5)
sfc_levels = np.arange(0, 3.10, 0.10)
swh_levels = np.arange(4, 14.5, 0.5)
sst_levels = np.arange(28.0, 29.65, 0.05)
stress_levels = np.arange(0,5.2,0.2)
stress_diff_levels = np.arange(-1.0,1.1,0.1)

#Vector Settings
skip = 15

#Labels
ctl_labels = ['(a)', '(d)', '(g)', '(j)']
exp_labels = ['(b)', '(e)', '(h)', '(k)']
diff_labels = ['(c)', '(f)', '(i)', '(l)']


time_indices = np.arange(1, qual_data.sizes['time'], 2)  # use step=1 for every time
nrows = len(time_indices)

clevs, clist = anom_cbar(20, len(swh_diff_levels), 'blue', 'tomato')


fig, axes = plt.subplots(nrows=nrows, ncols=3, figsize=(12.5, 2.6 * nrows))
if nrows == 1:
    axes = np.array([axes])

for r, i in enumerate(time_indices):
    print('Working on ' + str(qual_data.time.isel(time=i).values))

    

    ax1 = axes[r, 0]
    ax2 = axes[r, 1]
    ax3 = axes[r, 2]

    awo_wind = ax1.contourf(qual_data.x, qual_data.y, qual_data.swh_awo_qc.isel(time=i), 
                            levels=swh_levels, cmap=swh_cmaps, extend='both')
    # awo_wind.cmap.set_bad('saddlebrown')
    mask_ctl = np.isnan(qual_data.swh_awo_qc.isel(time=i))
    if np.any(mask_ctl):
        mask_ma = np.ma.masked_where(~mask_ctl, mask_ctl)  # mask only True locations
        ax1.pcolormesh(qual_data.x, qual_data.y, mask_ma, cmap=matplotlib.colors.ListedColormap(['saddlebrown']),
                        shading='auto', zorder=-10)

    
    awo_ws_wind = ax2.contourf(qual_data.x, qual_data.y, qual_data.swh_awo_ws_qc.isel(time=i),
                              levels=swh_levels, cmap=swh_cmaps, extend='both')
    mask_exp = np.isnan(qual_data.swh_awo_ws_qc.isel(time=i))
    if np.any(mask_exp):
        mask_ma = np.ma.masked_where(~mask_exp, mask_exp)  # mask only True locations
        ax2.pcolormesh(qual_data.x, qual_data.y, mask_ma, cmap=matplotlib.colors.ListedColormap(['saddlebrown']),
                        shading='auto', zorder=-10)
    
    swh_diff = ax3.contourf(qual_data.x, qual_data.y, qual_data.swh_awo_qc.isel(time=i) - qual_data.swh_awo_ws_qc.isel(time=i),
                             levels=swh_diff_levels, colors=clist, extend='both')
    mask_diff = np.isnan(qual_data.swh_awo_qc.isel(time=i) - qual_data.swh_awo_ws_qc.isel(time=i))
    if np.any(mask_diff):
        mask_ma = np.ma.masked_where(~mask_diff, mask_diff)  # mask only True locations
        ax3.pcolormesh(qual_data.x, qual_data.y, mask_ma, cmap=matplotlib.colors.ListedColormap(['saddlebrown']),
                        shading='auto', zorder=-10)
    
    add_storm_relative_essentials(ax1)
    add_storm_relative_axis_labels(ax1, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
    add_corner_label(ax1, x_pos, y_pos, ctl_labels[r], fontsize=fontsize)
    ax1.set_aspect('equal')
    ax1.set_title(f'$CTL$ {mid_periods[i]:%Y-%m-%d %H}', fontsize=fontsize, pad=1)

    add_storm_relative_essentials(ax2)
    add_storm_relative_axis_labels(ax2, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
    add_corner_label(ax2, x_pos, y_pos, exp_labels[r], fontsize=fontsize)
    ax2.set_aspect('equal')
    ax2.set_title(f'$EXP$ {mid_periods[i]:%Y-%m-%d %H}', fontsize=fontsize, pad=1)


    add_storm_relative_essentials(ax3)
    add_storm_relative_axis_labels(ax3, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
    add_corner_label(ax3, x_pos, y_pos, diff_labels[r], fontsize=fontsize)
    ax3.set_aspect('equal')
    ax3.set_title(f'$CTL - EXP$ {mid_periods[i]:%Y-%m-%d %H}', fontsize=fontsize, pad=1)

fig.subplots_adjust(right=0.86)
cbar_diff = fig.add_axes([0.98, 0.05, 0.02, 0.90])
cb_diff = fig.colorbar(swh_diff, cax=cbar_diff, ticks=swh_diff_levels, shrink=1, aspect=aspect)
cb_diff.ax.tick_params(color='k', length=1, width=0.5, labelsize=12, pad=0.002)
cb_diff.set_label('(CTL - EXP) Significant Wave Height (m)', fontsize=12, labelpad=2)

cbar_ax = fig.add_axes([0.88, 0.05, 0.02, 0.90])
cb = fig.colorbar(awo_wind, cax=cbar_ax, ticks=swh_levels, shrink=1, aspect=aspect)
cb.ax.tick_params(color='k', length=1, width=0.5, labelsize=12, pad=0.002)
cb.set_label('Significant Wave Height (m)', fontsize=12, labelpad=2)


fig.tight_layout(pad=0.1, w_pad=-15, h_pad=1)

fig.savefig(PNG + 'fig9_significant_wave_height_comparison_storm_radii_ESS.png', dpi=300, bbox_inches='tight')