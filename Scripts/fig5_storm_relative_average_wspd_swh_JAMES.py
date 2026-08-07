import xarray as xr
import cmaps
import numpy as np
import matplotlib.pyplot as plt
import os

from helpers import *

def plot_storm_radius_34(ax, row, prefix='CTL', color='k', lw=2, ls='-'):
        # Plot a smooth contour by concatenating all quadrants for the given prefix (CTL or EXP)
        theta_full = np.concatenate([
            np.linspace(0, np.pi/2, 100),                # NE
            np.linspace(np.pi/2, np.pi, 100),            # NW
            np.linspace(np.pi, 3*np.pi/2, 100),          # SW
            np.linspace(3*np.pi/2, 2*np.pi, 100)         # SE
        ])
        radius_full = np.concatenate([
            np.full(100, row[f'{prefix}_R34_NE']),
            np.full(100, row[f'{prefix}_R34_NW']),
            np.full(100, row[f'{prefix}_R34_SW']),
            np.full(100, row[f'{prefix}_R34_SE'])
        ])
        x_full = radius_full * np.cos(theta_full)
        y_full = radius_full * np.sin(theta_full)
        ax.plot(x_full, y_full, color=color, lw=lw, ls=ls, zorder=10)

def plot_storm_radius_50(ax, row, prefix='CTL', color='k', lw=2, ls='-'):
        # Plot a smooth contour by concatenating all quadrants for the given prefix (CTL or EXP)
        theta_full = np.concatenate([
            np.linspace(0, np.pi/2, 100),                # NE
            np.linspace(np.pi/2, np.pi, 100),            # NW
            np.linspace(np.pi, 3*np.pi/2, 100),          # SW
            np.linspace(3*np.pi/2, 2*np.pi, 100)         # SE
        ])
        radius_full = np.concatenate([
            np.full(100, row[f'{prefix}_R50_NE']),
            np.full(100, row[f'{prefix}_R50_NW']),
            np.full(100, row[f'{prefix}_R50_SW']),
            np.full(100, row[f'{prefix}_R50_SE'])
        ])
        x_full = radius_full * np.cos(theta_full)
        y_full = radius_full * np.sin(theta_full)
        ax.plot(x_full, y_full, color=color, lw=lw, ls=ls, zorder=10)

def plot_storm_radius_64(ax, row, prefix='CTL', color='k', lw=2, ls='-'):
        # Plot a smooth contour by concatenating all quadrants for the given prefix (CTL or EXP)
        theta_full = np.concatenate([
            np.linspace(0, np.pi/2, 100),                # NE
            np.linspace(np.pi/2, np.pi, 100),            # NW
            np.linspace(np.pi, 3*np.pi/2, 100),          # SW
            np.linspace(3*np.pi/2, 2*np.pi, 100)         # SE
        ])
        radius_full = np.concatenate([
            np.full(100, row[f'{prefix}_R64_NE']),
            np.full(100, row[f'{prefix}_R64_NW']),
            np.full(100, row[f'{prefix}_R64_SW']),
            np.full(100, row[f'{prefix}_R64_SE'])
        ])
        x_full = radius_full * np.cos(theta_full)
        y_full = radius_full * np.sin(theta_full)
        ax.plot(x_full, y_full, color=color, lw=lw, ls=ls, zorder=10)

# Data_Dir = '/home/disk/orca/adaley17/Research/Stress_Separation/Hurricane_Earl/Data/'
DATA_PATH = '/home/disk/orca/adaley17/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/notebooks/'
PATH = '/home/disk/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/'
PNG='/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'
file = 'storm_relative_qualit_controlled_data.2010090106_2010090111.nc'

data = xr.open_mfdataset(PATH + file)

awo_track_file = 'earl_5.log_awo.csv'
awo_ws_track_file = 'earl_5.log_awo-ws.csv'
radius_file = 'radius_df.csv'

awo_csv_file_path = os.path.join(PATH, awo_track_file)
awo_ws_csv_file_path = os.path.join(PATH, awo_ws_track_file)

#Radius DataFrame
radius_df = pd.read_csv(DATA_PATH + radius_file)

# Read the CSV file into a DataFrame
awo_track_df = pd.read_csv(awo_csv_file_path) #AWO
awo_ws_track_df = pd.read_csv(awo_ws_csv_file_path) #AWO_ws

#Compute Wind Radii
start_time = '2010-09-01 06:00:00' #Enter Start Time Here
end_time = '2010-09-01 11:00:00' #Enter End Time Here
time_of_int = '2010-09-01 09:00:00' #Enter Time of Interest Here
# print(awo_track_df.head())
# subset_awo_track_df, max_r34_across_time, max_r50_across_time, max_r64_across_time = compute_wind_radii(awo_track_df, start_time, end_time)
# subset_awo_ws_track_df, max_r34_across_time_ws, max_r50_across_time_ws, max_r64_across_time_ws = compute_wind_radii(awo_ws_track_file, start_time, end_time)
awo_max_r64, awo_max_r50, awo_max_r34 = compute_max_wind_radii(awo_track_df, time_of_int)
awo_ws_r64, awo_ws_r50, awo_ws_r34 = compute_max_wind_radii(awo_ws_track_df, time_of_int)

# print(subset_awo_track_df[['r34_ne', 'r34_se', 'r34sw', 'r34nw']])

print('R64', awo_max_r64, awo_ws_r64)
print('R50', awo_max_r50, awo_ws_r50)
print('R34', awo_max_r34, awo_ws_r34)

# # Compute the average radius for r34, r50, and r64
# subset_awo_track_df['r34_max'] = subset_awo_track_df[['r34_ne', 'r34_se', 'r34sw', 'r34nw']].max(axis=1)
# subset_awo_track_df['r50_max'] = subset_awo_track_df[['r50ne', 'r50se', 'r50sw', 'r50nw']].max(axis=1)
# subset_awo_track_df['r64_max'] = subset_awo_track_df[['r64ne', 'r64se', 'r64sw', 'r64nw']].max(axis=1)

# # Compute the average wind radii across time
# max_r34_across_time = subset_awo_track_df['r34_max'].values.max() 
# max_r50_across_time = subset_awo_track_df['r50_max'].values.max()
# max_r64_across_time = subset_awo_track_df['r64_max'].values.max()

# print(subset_awo_track_df['r34_max'])

theta = np.arctan2(data['y'], data['x'])
R = np.sqrt(data['x']**2 + data['y']**2)


#Plotting Options
xmin = -250
xmax = 250
ymin = xmin
ymax = xmax
fontsize=4.5
labelsize=4
labelpad=0.25
zorder=-5
zorder_new=-1
y_pos = 0.87
x_pos = 0.085
zorder=50

ticks=np.arange(-250,300,50)

#Color Maps
wspd_cmaps = cmaps.cmp_haxby_r
swh_cmap = cmaps.NMCRef
diff_cmap = cmaps.MPL_bwr
sst_cmap = cmaps.MPL_jet
sfc_cmap=cmaps.BkBlAqGrYeOrReViWh200_r

#Levels
wspd_diff_levels = np.arange(-5,6,1)
swh_diff_levels = np.arange(-0.7,0.8,0.1) 
sst_diff_levels = np.arange(-0.40,0.45, 0.05)
sfc_diff_levels = np.arange(-0.5,0.55,0.05)
wspd_levels = np.arange(7,55,2.5)
sfc_levels = np.arange(0, 3.10, 0.10)
swh_levels = np.arange(4, 14.5, 0.5)
sst_levels = np.arange(28.0, 29.65, 0.05)

clevs_wspd, clist_wspd = anom_cbar(20, len(wspd_diff_levels), 'blue', 'tomato')
clevs_swh, clist_swh = anom_cbar(0.5, len(swh_diff_levels), 'blue', 'tomato')

#Colorbar Options
aspect = 12
shrink = 0.30

gridsize = (2, 3)
fig = plt.figure(figsize=(4, 4))

x_reshape, y_reshpae = np.meshgrid(data['x'], data['y'])

# AWO 
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)

awo_wspd = ax1.contourf(data['x'], data['y'], data['wspd_awo_qc'][0], 
                        levels=wspd_levels, cmap=wspd_cmaps, extend='both')

# CS = ax1.contour(data['x'], data['y'],data['wspd_awo_ws_qc'][0], [17], linewidths=0.5, colors='white')
# ax1.clabel(CS, inline=True, fontsize=5)

#Wind Radii
# CS_34 = plt.contour(R * np.cos(theta), R * np.sin(theta), R, levels=[awo_max_r34], colors='white', zorder=zorder, linewidths=1)
# CS_50 = plt.contour(R * np.cos(theta), R * np.sin(theta), R, levels=[awo_max_r50], colors='magenta', zorder=zorder, linewidths=1)
# CS_64 = plt.contour(R * np.cos(theta), R * np.sin(theta), R, levels=[awo_max_r64], colors='blue', zorder=zorder, linewidths=1)
# plt.clabel(CS_34, inline=True, fontsize=15, fmt='%1.0f km^2')
# plt.clabel(CS_50, inline=True, fontsize=15, fmt='%1.0f km^2')
# plt.clabel(CS_64, inline=True, fontsize=15, fmt='%1.0f km^2')

# ax1.text(-50, -115, 'R64 = {:.2f} $km^2$'.format(max_r64_across_time), color='blue', fontsize=fontsize)
# ax1.text(-50, -140, 'R50 = {:.2f} $km^2$'.format(max_r50_across_time), color='magenta', fontsize=fontsize)
# ax1.text(-50, -165, 'R34 = {:.2f} $km^2$'.format(max_r34_across_time), color='white', fontsize=fontsize)


hcb = fig.colorbar(awo_wspd, shrink=shrink, aspect=aspect, ax=ax1, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=labelsize, pad=0.002)

add_storm_relative_essentials(ax1)
add_storm_relative_axis_labels(ax1, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax1, x_pos, y_pos, '(a)', fontsize=fontsize)
ax1.set_title('$CTL$ $U_{10}$ ($m/s$)', fontsize=fontsize, pad=1)
ax1.set_aspect('equal')

#AWO_ws 
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)

awo_ws_wspd = ax2.contourf(data['x'], data['y'],data['wspd_awo_ws_qc'][0], 
                        levels=wspd_levels, cmap=wspd_cmaps, extend='both')

# CS = ax2.contour(data['x'], data['y'],data['wspd_awo_ws_qc'][0], [17], linewidths=0.5, colors='white')
# ax2.clabel(CS, inline=True, fontsize=5)

# CS_34_ws = plt.contour(R * np.cos(theta), R * np.sin(theta), R, levels=[awo_ws_r34], colors='white', zorder=zorder, linewidths=1)
# CS_50_ws = plt.contour(R * np.cos(theta), R * np.sin(theta), R, levels=[awo_ws_r50], colors='magenta', zorder=zorder, linewidths=1)
# CS_64_ws = plt.contour(R * np.cos(theta), R * np.sin(theta), R, levels=[awo_max_r64], colors='blue', zorder=zorder, linewidths=1)
# plt.clabel(CS_34_ws, inline=True, fontsize=15, fmt='%1.0f $km^2$')
# plt.clabel(CS_50_ws, inline=True, fontsize=15, fmt='%1.0f $km^2$')
# plt.clabel(CS_64_ws, inline=True, fontsize=15, fmt='%1.0f $km^2$')

# ax2.text(-50, -115, 'R64 = {:.2f} $km^2$'.format(max_r64_across_time_ws), color='blue', fontsize=fontsize)
# ax2.text(-50, -140, 'R50 = {:.2f} $km^2$'.format(max_r50_across_time_ws), color='magenta', fontsize=fontsize)
# ax2.text(-50, -165, 'R34 = {:.2f} $km^2$'.format(max_r34_across_time_ws), color='white', fontsize=fontsize)

hcb = fig.colorbar(awo_ws_wspd, shrink=shrink, aspect=aspect, ax=ax2, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=labelsize, pad=0.002)

add_storm_relative_essentials(ax2)
# add_axis_labels(axis, fontsize, labelpad, labelsize, xticks, yticks):
add_storm_relative_axis_labels(ax2, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax2, x_pos, y_pos, '(b)', fontsize=fontsize)
ax2.set_title('$EXP$ $U_{10}$ ($m/s$)', fontsize=fontsize, pad=1)
ax2.set_aspect('equal')

# Difference
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)


wspd_awo_awo_ws = ax3.contourf(data['x'], data['y'], data['wspd_awo_qc'][0] - data['wspd_awo_ws_qc'][0], 
                        levels=wspd_diff_levels, colors=clist_wspd, extend='both')

row = radius_df[radius_df['time'] == time_of_int]
# Plot CTL line (black, solid)
plot_storm_radius_34(ax3, row, prefix='CTL', color='k', lw=1, ls='-')
plot_storm_radius_34(ax3, row, prefix='EXP', color='k', lw=1, ls='--')

plot_storm_radius_50(ax3, row, prefix='CTL', color='orange', lw=1, ls='-')
plot_storm_radius_50(ax3, row, prefix='EXP', color='orange', lw=1, ls='--')

plot_storm_radius_64(ax3, row, prefix='CTL', color='green', lw=1, ls='-')
plot_storm_radius_64(ax3, row, prefix='EXP', color='green', lw=1, ls='--')


hcb = fig.colorbar(wspd_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax3, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=labelsize, pad=0.002)

add_storm_relative_essentials(ax3)
add_storm_relative_axis_labels(ax3, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax3, x_pos, y_pos, '(c)', fontsize=fontsize)
ax3.set_title('$CTL$ - $EXP$ $U_{10}$ ($m/s$)', fontsize=fontsize, pad=1)
ax3.set_aspect('equal')

#SWH
#AWO
ax4 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)

awo_swh = ax4.contourf(data['x'], data['y'], data['swh_awo_qc'][0], 
                        levels=swh_levels, cmap=swh_cmap, extend='both')
# depth = ax4.contour(awo_x_dist_6, awo_y_dist_6, bathy,  linewidths=0.5, levels=depths, colors='black')

hcb = fig.colorbar(awo_swh, shrink=shrink, aspect=aspect, ax=ax4, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax4)
add_storm_relative_axis_labels(ax4, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax4, x_pos, y_pos, '(d)', fontsize=fontsize)

ax4.set_title('$CTL$ $H_{s}$ ($m$)', fontsize=fontsize, pad=1)
ax4.set_aspect('equal')

#AWO_ws
ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)

awo_ws_swh = ax5.contourf(data['x'], data['y'], data['swh_awo_ws_qc'][0], 
                        levels=swh_levels, cmap=swh_cmap, extend='both')

hcb = fig.colorbar(awo_ws_swh, shrink=shrink, aspect=aspect, ax=ax5, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax5)
add_storm_relative_axis_labels(ax5, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax5, x_pos, y_pos, '(e)', fontsize=fontsize)
ax5.set_title('$EXP$ $H_{s}$ ($m$)', fontsize=fontsize, pad=1)
ax5.set_aspect('equal')



#Diff
ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)

swh_awo_awo_ws = ax6.contourf(data['x'], data['y'], data['swh_awo_qc'][0] - data['swh_awo_ws_qc'][0], 
                        levels=swh_diff_levels, colors=clist_swh, extend='both')

hcb = fig.colorbar(swh_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax6, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax6)
add_storm_relative_axis_labels(ax6, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax6, x_pos, y_pos, '(f)', fontsize=fontsize)
ax6.set_title('$CTL$ - $EXP$ $H_{s}$ ($m$)', fontsize=fontsize, pad=1)
ax6.set_aspect('equal')


fig.tight_layout(pad=0.01, w_pad=0.5, h_pad=-9.5)

plt.savefig(PNG + 'storm_relative_average_wspd_swh_wind_radii_JAMES.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)
