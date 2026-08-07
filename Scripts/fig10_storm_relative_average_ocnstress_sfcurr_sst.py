import xarray as xr
import cmaps
import numpy as np
import matplotlib.pyplot as plt
import os

from helpers import *

Data_Dir = '/home/disk/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/notebooks/'
PATH = '/home/disk/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/'
PNG='/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'
file = 'storm_relative_qualit_controlled_data.2010090106_2010090111.nc'

data = xr.open_mfdataset(Data_Dir + file)

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

ticks=np.arange(-250,300,50)

#Color Maps
wspd_cmaps = cmaps.MPL_gist_rainbow_r
swh_cmap = cmaps.NMCRef
diff_cmap = cmaps.MPL_bwr
sst_cmap = cmaps.MPL_jet
sfc_cmap = cmaps.WhiteBlueGreenYellowRed
stress_cmap = cmaps.precip3_16lev

#Levels
wspd_diff_levels = np.arange(-5,6,1)
swh_diff_levels = np.arange(-1,1.1,0.1) 
sst_diff_levels = np.arange(-0.40,0.45, 0.05)
sfc_diff_levels = np.arange(-0.5,0.55,0.05)
wspd_levels = np.arange(5,47.5,2.5)
sfc_levels = np.arange(0.4, 2.05, 0.05)
swh_levels = np.arange(4, 14.5, 0.5)
sst_levels = np.arange(28.0, 29.65, 0.05)
stress_levels = np.arange(0,5.2,0.2)
stress_diff_levels = np.arange(-1.0,1.1,0.1)

clevs_sst, clist_sst = anom_cbar(20, len(sst_diff_levels), 'blue', 'tomato')
clevs_stress, clist_stress = anom_cbar(20, len(stress_diff_levels), 'blue', 'tomato')
clevs_sfc, clist_sfc = anom_cbar(20, len(sfc_diff_levels), 'blue', 'tomato')

#Colorbar Options
aspect = 12
shrink = 0.50

gridsize = (3, 3)
fig = plt.figure(figsize=(4, 4))

x_reshape, y_reshpae = np.meshgrid(data['x'], data['y'])

#Ocean Stress
#AWO
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)

awo_stress = ax1.contourf(data['x'], data['y'], data['awo_ocn_stress'], 
                        cmap=stress_cmap,levels=stress_levels, extend='max')

hcb = fig.colorbar(awo_stress, shrink=shrink, aspect=aspect, ax=ax1, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax1)
add_storm_relative_axis_labels(ax1, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax1, x_pos, y_pos, '(a)', fontsize=fontsize)
ax1.set_title('$CTL$ $\\tau_{ocn}$ ($N/m^{2}$)', fontsize=fontsize, pad=1)
ax1.set_aspect('equal')

#AWO_ws
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
awo_ws_stress = ax2.contourf(data['x'], data['y'], data['awo_ws_ocn_stress'], 
                        cmap=stress_cmap,levels=stress_levels, extend='max')

hcb = fig.colorbar(awo_ws_stress, shrink=shrink, aspect=aspect, ax=ax2, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax2)
add_storm_relative_axis_labels(ax2, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax2, x_pos, y_pos, '(b)', fontsize=fontsize)
ax2.set_title('$EXP$ $\\tau_{ocn}$ ($N/m^{2}$)', fontsize=fontsize, pad=1)
ax2.set_aspect('equal')

#Diff
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)  
stress_diff = ax3.contourf(data['x'], data['y'], data['awo_ocn_stress'] - data['awo_ws_ocn_stress'],
                        colors=clist_stress, levels=stress_diff_levels, extend='both')

hcb = fig.colorbar(stress_diff, shrink=shrink, aspect=aspect, ax=ax3, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax3)
add_storm_relative_axis_labels(ax3, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax3, x_pos, y_pos, '(c)', fontsize=fontsize)
ax3.set_title('$CTL$ - $EXP$ $\\tau_{ocn}$ ($N/m^{2}$)', fontsize=fontsize, pad=1)
ax3.set_aspect('equal')

#Surface Currents
#AWO
ax4 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)

awo_sfc = ax4.contourf(data['x'], data['y'], data['sfc_awo_qc'][0],
                        cmap=sfc_cmap,levels=sfc_levels, extend='both')

hcb = fig.colorbar(awo_sfc, shrink=shrink, aspect=aspect, ax=ax4, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax4)
add_storm_relative_axis_labels(ax4, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax4, x_pos, y_pos, '(d)', fontsize=fontsize)
ax4.set_title('$CTL$ $SFC$ ($m/s$)', fontsize=fontsize, pad=1)
ax4.set_aspect('equal')

#AWO_ws
ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)

awo_ws_sfc = ax5.contourf(data['x'], data['y'], data['sfc_awo_ws_qc'][0], 
                        cmap=sfc_cmap,levels=sfc_levels, extend='both')


hcb = fig.colorbar(awo_ws_sfc, shrink=shrink, aspect=aspect, ax=ax5, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax5)
add_storm_relative_axis_labels(ax5, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax5, x_pos, y_pos, '(e)', fontsize=fontsize)
ax5.set_title('$EXP$ $SFC$ ($m/s$)', fontsize=fontsize, pad=1)
ax5.set_aspect('equal')

#Diff
ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)

sfc_awo_awo_ws = ax6.contourf(data['x'], data['y'], data['sfc_awo_qc'][0] - data['sfc_awo_ws_qc'][0], 
                        colors=clist_sfc, levels=sfc_diff_levels, extend='both')

CS = ax6.contour(x_reshape, y_reshpae, data['sst_awo_ws_qc'][0], [28.5], linewidths=0.5, colors='magenta')
# ax9.axhline(y=-98, xmin=-50, xmax=50, color='yellow', linestyle='-')
ax6.plot((-100, 100), (-100, -100), color='yellow', linestyle='-') 

hcb = fig.colorbar(sfc_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax6, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax6)
add_storm_relative_axis_labels(ax6, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax6, x_pos, y_pos, '(f)', fontsize=fontsize)
ax6.set_title('$CTL$ - $EXP$ $SFC$ ($m/s$)', fontsize=fontsize, pad=1)
ax6.set_aspect('equal')

#Surface Currents
#AWO
ax7 = plt.subplot2grid(gridsize, (2, 0), colspan=1, rowspan=1)

awo_sst = ax7.contourf(data['x'], data['y'], data['sst_awo_qc'][0], 
                         levels=sst_levels, cmap=sst_cmap, extend='both', zorder=zorder)
CS = ax7.contour(x_reshape, y_reshpae, data['sst_awo_ws_qc'][0], [28.5], linewidths=0.5, colors='magenta')


hcb = fig.colorbar(awo_sst, shrink=shrink, aspect=aspect, ax=ax7, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax7)
add_storm_relative_axis_labels(ax7, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax7, x_pos, y_pos, '(g)', fontsize=fontsize)
ax7.set_title('$CTL$ $T$ ($^{\circ}C$)', fontsize=fontsize, pad=1)
ax7.set_aspect('equal')


#AWO_ws
ax8 = plt.subplot2grid(gridsize, (2, 1), colspan=1, rowspan=1)

awo_ws_sst = ax8.contourf(data['x'], data['y'], data['sst_awo_ws_qc'][0], 
                        levels=sst_levels, cmap=sst_cmap, extend='both')
CS = ax8.contour(x_reshape, y_reshpae, data['sst_awo_ws_qc'][0], [28.5], linewidths=0.5, colors='magenta')

hcb = fig.colorbar(awo_ws_sst, shrink=shrink, aspect=aspect, ax=ax8, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax8)
add_storm_relative_axis_labels(ax8, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax8, x_pos, y_pos, '(h)', fontsize=fontsize)
ax8.set_title('$EXP$ $T$ ($^{\circ}C$)', fontsize=fontsize, pad=1)
ax8.set_aspect('equal')

# Diff
ax9 = plt.subplot2grid(gridsize, (2, 2), colspan=1, rowspan=1)

sst_awo_awo_ws = ax9.contourf(data['x'], data['y'], data['sst_awo_qc'][0] - data['sst_awo_ws_qc'][0], 
                         colors=clist_sst, levels=sst_diff_levels, extend='both', zorder=zorder)

# print(np.nanmax(data['sst_awo_qc'][0] - data['wspd_awo_ws_qc'][0]))

# CS = ax12.contour(x_reshape, y_reshpae, data['sst_awo_qc'][0] - data['sst_awo_ws_qc'][0], [0.20], linewidths=0.5, colors='yellow')
CS = ax9.contour(x_reshape, y_reshpae, data['sst_awo_ws_qc'][0], [28.5], linewidths=0.5, colors='magenta')
ax9.plot((-100, 100), (-100, -100), color='yellow', linestyle='-')

hcb = fig.colorbar(sst_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax9, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax9)
add_storm_relative_axis_labels(ax9, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax9, x_pos, y_pos, '(i)', fontsize=fontsize)
ax9.set_title('$CTL$ - $EXP$ $T$ ($^{\circ}C$)', fontsize=fontsize, pad=1)
ax9.set_aspect('equal')

fig.tight_layout(pad=0.01, w_pad=0.5, h_pad=-4)

plt.savefig(PNG + 'fig10_storm_relative_average_ocnstress_sfcurr_sst_ESS.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)
