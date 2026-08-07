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
from scipy import interpolate
# from mompy import *
from helpers import *
from matplotlib.patches import Wedge



DATA_PATH = '/home/disk/orca/adaley17/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/notebooks/'
SUM_PATH = '/home/disk/orca/adaley17/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/'
PNG = '/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/' # Path to save figures
data = 'storm_relative_qualit_controlled_data.*.nc'

radius_df = pd.read_csv(DATA_PATH + 'radius_df.csv', parse_dates=['time'])

# print(radius_df)

qual_data = xr.open_mfdataset(DATA_PATH + data, combine='by_coords')
ctl_mod_radius_data = extract_storm_radius_data(SUM_PATH + 'earl_5.log_awo.csv')
exp_mod_radius_data = extract_storm_radius_data(SUM_PATH + 'earl_5.log_awo-ws.csv')


time_indices = np.arange(1, qual_data.sizes['time'], 2)  # use step=1 for every time

mid_periods = pd.date_range(start='2010-08-30T15:00:00', end='2010-09-01T09:00:00',
                            freq='6h') # Define mid-periods for analysis

radius_r50_df = radius_df[['time'] + [col for col in radius_df.columns if 'R50' in col]]
r50_cols = [col for col in radius_r50_df.columns if col != 'time']

for _, row in radius_r50_df.iterrows():
    print(f"{row['time']:%Y-%m-%d %H:%M} -> max R50: {row[r50_cols].max():.3f}")
    if _ == 0:
        max_r50_rows = []

    max_col = row[r50_cols].idxmax()
    max_r50_rows.append({
        'time': row['time'],
        'max_R50': row[r50_cols].max(),
        'max_R50_source': max_col
    })

    radius_r50_max_df = pd.DataFrame(max_r50_rows)
    
    print(radius_r50_max_df)

for r, i in enumerate(time_indices):
    print('Working on ' + str(qual_data.time.isel(time=i).values))
    # ax = axes[r]
    # Create a mask for points within the storm's R50 (50-kt wind) radius
    xx, yy = np.meshgrid(qual_data.x, qual_data.y)
    mask = np.sqrt(xx**2 + yy**2) < radius_r50_max_df.loc[radius_r50_max_df['time'] == mid_periods[i], 'max_R50'].values[0]
    diff = qual_data.wspd_awo_qc.isel(time=i) - qual_data.wspd_awo_ws_qc.isel(time=i)
    ctl = qual_data.wspd_awo_qc.isel(time=i).where(~mask)
    exp = qual_data.wspd_awo_ws_qc.isel(time=i).where(~mask)
    diff = diff.where(~mask)

    # Collect results for each time step in a list
    if r == 0:
        results_list = []
    results_list.append({
        'Time': str(qual_data.time.isel(time=i).values),
        'CTL mean': ctl.mean().compute().item(),
        'EXP mean': exp.mean().compute().item(),
        'DIFF mean': diff.mean().compute().item(),
    })

# After the loop, convert to DataFrame
results_df = pd.DataFrame(results_list)

gridsize = (2, 1)

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(4, 4))

axes[0].scatter(mid_periods[1:][::2], results_df['CTL mean'], label='CTL', color='red')
axes[0].scatter(mid_periods[1:][::2], results_df['EXP mean'], label='EXP', color='cyan')


axes[0].set_ylabel('Mean Wind \nSpeed (m/s)')
axes[0].set_xticks(mid_periods[1:][::2], labels=[f'{t:%m-%d %H}' for t in mid_periods[1:][::2]], rotation=45)
axes[0].set_xticklabels([])
axes[0].grid(linestyle='--', alpha=0.5)

axes[1].scatter(mid_periods[1:][::2], results_df['DIFF mean'], label='DIFF', color='purple')
axes[1].set_xticks(mid_periods[1:][::2], labels=[f'{t:%m-%d %H}' for t in mid_periods[1:][::2]], rotation=45)
axes[1].hlines(0, mid_periods[1:][::2][0], mid_periods[1:][::2][-1], colors='gray', linestyles='-', alpha=1, lw=2.5)
axes[1].set_ylabel('Mean Wind Speed \nDifference (m/s)')
axes[1].grid(linestyle='--', alpha=0.5)

axes[1].set_xlabel('Time (mm-dd HH)')

axes[0].legend(loc='upper left', shadow=True)

plt.savefig(PNG + 'fig8_avg_50kt_wspd_diff_ctl_vs_exp_ESS.png', dpi=300, bbox_inches='tight')

# plt.show()