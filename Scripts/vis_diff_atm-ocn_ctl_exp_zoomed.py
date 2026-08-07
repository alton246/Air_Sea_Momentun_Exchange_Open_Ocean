import xarray as xr
import cmaps
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

from helpers import *





Data_Dir = '/home/disk/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/notebooks/'
PATH = '/home/disk/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/'
PNG='/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'
file = 'storm_relative_qualit_controlled_data.*.nc'

data = xr.open_mfdataset(Data_Dir + file)
data['ctl_stress_diff'] = data['awo_atm_stress'] - data['awo_ocn_stress']
data['exp_stress_diff'] = data['awo_ws_atm_stress'] - data['awo_ws_ocn_stress']

data['awo_land_mask'] = data['ctl_stress_diff'].where(~np.isnan(data['ctl_stress_diff']), 999)
data['awo_ws_land_mask'] = data['exp_stress_diff'].where(~np.isnan(data['exp_stress_diff']), 999)

# print(data['awo_land_mask'].values)


awo_track_file = 'earl_5.log_awo.csv'
awo_ws_track_file = 'earl_5.log_awo-ws.csv'

awo_csv_file_path = os.path.join(PATH, awo_track_file)
awo_ws_csv_file_path = os.path.join(PATH, awo_ws_track_file)

# Read the CSV file into a DataFrame
awo_track_df = pd.read_csv(awo_csv_file_path) #AWO
awo_ws_track_df = pd.read_csv(awo_ws_csv_file_path) #AWO_ws

#Compute Wind Radii
start_time = '2010-09-01 06:00:00' #Enter Start Time Here
end_time = '2010-09-01 11:00:00' #Enter End Time Here
time_of_int = '2010-09-01 09:00:00' #Enter Time of Interest Here

r = np.sqrt(data.x**2 + data.y**2) # Compute radial distance
dist = r.values.flatten() # Flatten distance array

awo_stress_ratio_bins = np.linspace(0, 1.6, 50) # Define stress ratio bins
dist_bins = np.linspace(20, 200, 100) # Define distance bins

awo_stress_ratio_mean_comb = [] # Initialize list to store mean CTL stress ratios
awo_ws_stress_ratio_mean_comb = [] # Initialize list to store EXP stress ratios


mid_periods = pd.date_range(start='2010-08-30T15:00:00', end='2010-09-01T09:00:00', 
                            freq='6h') # Define mid-periods for analysis
n = len(mid_periods) # Number of mid-periods

#Plotting Options
xmin = -250
xmax = 250
ymin = xmin
ymax = xmax
fontsize=4
labelsize=3.75
labelpad=0.25
zorder=-5
zorder_new=-1
y_pos = 0.87
x_pos = 0.085
lw=0.8
alpha=0.5
linewidth=0.5

colors=['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet', 'k']
xlim = (0, 200)
ylim = (0.80, 1.05)
zorder=50

ticks=np.arange(-250,300,50)
ticks_zoomed = np.arange(-50,60,10)


#Color Maps
wspd_cmaps = cmaps.MPL_gist_rainbow_r
swh_cmap = cmaps.NMCRef
diff_cmap = cmaps.BlWhRe
sst_cmap = cmaps.MPL_jet
sfc_cmap=cmaps.BkBlAqGrYeOrReViWh200_r
stress_cmap = cmaps.GMT_no_green

#Levels
wspd_diff_levels = np.arange(-5,6,1)
swh_diff_levels = np.arange(-1,1.1,0.1) 
# sst_diff_levels = np.arange(-0.20,0.225, 0.025)
sst_diff_levels = np.arange(-0.10,0.125, 0.025)
sfc_diff_levels = np.arange(-0.5,0.55,0.05)
wspd_levels = np.arange(5,47.5,2.5)
sfc_levels = np.arange(0, 3.10, 0.10)
swh_levels = np.arange(4, 14.5, 0.5)
sst_levels = np.arange(28.0, 29.65, 0.05)
stress_levels = np.arange(0.75,1.30,0.05)
stress_diff_levels = np.arange(-1.0,1.1,0.1)

# print(stress_levels)

clevs, clist = anom_cbar(20, len(sst_diff_levels), 'blue', 'tomato')

#Colorbar Options
aspect = 12
shrink = 0.98

gridsize = (4, 2)
fig = plt.figure(figsize=(4, 4))

x_reshape, y_reshpae = np.meshgrid(data['x'], data['y'])

print(data['time'].values[4])

#Ocean Stress
#AWO
axes = []
titles = []
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']
factors = [
    ('ctl_stress_diff', 1, '$CTL$'),
    ('ctl_stress_diff', 1, '$CTL$'),
    ('ctl_stress_diff', 3, '$CTL$'),
    ('ctl_stress_diff', 3, '$CTL$'),
    ('ctl_stress_diff', 5, '$CTL$'),
    ('ctl_stress_diff', 5, '$CTL$'),
    ('ctl_stress_diff', 7, '$CTL$'),
    ('ctl_stress_diff', 7, '$CTL$'),
]
positions = [(i, j) for i in range(4) for j in range(2)]

for idx, ((var, tidx, exp_label), pos) in enumerate(zip(factors, positions)):
    ax = plt.subplot2grid(gridsize, pos, colspan=1, rowspan=1)
    cf = ax.contourf(data['x'], data['y'], data[var][tidx], levels=sst_diff_levels, 
                     colors=clist, extend='both')
    
    cf.cmap.set_bad('saddlebrown')
    mask = np.isnan(data[var][tidx].values)
    if np.any(mask):
        mask_ma = np.ma.masked_where(~mask, mask)  # mask only True locations
        ax.pcolormesh(data['x'], data['y'], mask_ma, cmap=matplotlib.colors.ListedColormap(['saddlebrown']),
                        shading='auto', zorder=zorder_new)

    add_storm_relative_essentials(ax)
    add_storm_relative_axis_labels(ax, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
    add_corner_label(ax, x_pos, y_pos, labels[idx], fontsize=fontsize)
    ax.set_title(f'{exp_label}  {mid_periods[tidx].strftime("%Y-%m-%d %H")}',
                 fontsize=4.5, pad=1)
    ax.set_aspect('equal')
    axes.append(ax)

    if pos[1] == 1:
        ax.set_xlim(-50, 50)
        ax.set_ylim(-50, 50)
        ax.set_xticks(ticks_zoomed[::2])
        ax.set_yticks(ticks_zoomed[::2])
        ax.arrow(0, 0, 0, 15, width=0.5, color='k')

    axes.append(ax)
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.78, 0.05, 0.02, 0.90])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=sst_diff_levels,shrink=1, aspect=aspect)
cb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)
cb.set_label('ATM - OCN Stress Difference (N/m²)', fontsize=6, labelpad=2)


fig.tight_layout(pad=0, w_pad=-10, h_pad=0.5)

plt.savefig(PNG + 'atm-ocean_stress_diff_zoomed.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)
