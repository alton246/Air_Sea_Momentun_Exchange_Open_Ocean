import xarray as xr
import cmaps
import numpy as np
import matplotlib.pyplot as plt

from helpers import *

PATH = '/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentum_Exchange_TC_Coast/Notebooks/Figure2/'
PNG='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'
file = 'storm_relative_data.2010090100-2010090105.nc'

data = xr.open_mfdataset(PATH + file)


xmin = -250
xmax = 250
ymin = xmin
ymax = xmax
fontsize=3.75
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
sfc_cmap=cmaps.BkBlAqGrYeOrReViWh200_r

#Levels
wspd_diff_levels = np.arange(-5,6,1)
swh_diff_levels = np.arange(-1,1.1,0.1) 
sst_diff_levels = np.arange(-0.40,0.45, 0.05)
sfc_diff_levels = np.arange(-0.5,0.55,0.05)
wspd_levels = np.arange(5,47.5,2.5)
sfc_levels = np.arange(0, 3.10, 0.10)
swh_levels = np.arange(4, 14.5, 0.5)
sst_levels = np.arange(28.0, 29.65, 0.05)

#Colorbar Options
aspect = 12
shrink = 0.95

gridsize = (4, 3)
fig = plt.figure(figsize=(4, 4))


# AWO 
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)

awo_wspd = ax1.contourf(data['x'], data['y'], data['wspd_awo_qc'][0], 
                        levels=wspd_levels, cmap=wspd_cmaps, extend='both')

CS = ax1.contour(data['x'], data['y'],data['wspd_awo_ws_qc'][0], [17], linewidths=0.5, colors='white')
ax1.clabel(CS, inline=True, fontsize=5)


hcb = fig.colorbar(awo_wspd, shrink=shrink, aspect=aspect, ax=ax1, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=labelsize, pad=0.002)

add_storm_relative_essentials(ax1)
add_axis_labels(ax1, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax1, x_pos, y_pos, '(a)', fontsize=fontsize)
ax1.set_title('$AWO$-$CTL$ $U_{10}$ ($m/s$)', fontsize=fontsize, pad=1)
ax1.set_aspect('equal')

#AWO_ws 
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)

awo_ws_wspd = ax2.contourf(data['x'], data['y'],data['wspd_awo_ws_qc'][0], 
                        levels=wspd_levels, cmap=wspd_cmaps, extend='both')

CS = ax2.contour(data['x'], data['y'],data['wspd_awo_ws_qc'][0], [17], linewidths=0.5, colors='white')
ax2.clabel(CS, inline=True, fontsize=5)

hcb = fig.colorbar(awo_ws_wspd, shrink=shrink, aspect=aspect, ax=ax2, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=labelsize, pad=0.002)

add_storm_relative_essentials(ax2)
# add_axis_labels(axis, fontsize, labelpad, labelsize, xticks, yticks):
add_axis_labels(ax2, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax2, x_pos, y_pos, '(b)', fontsize=fontsize)
ax2.set_title('$AWO_{ws}$-$EXP$ $U_{10}$ ($m/s$)', fontsize=fontsize, pad=1)
ax2.set_aspect('equal')

# Difference
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)


wspd_awo_awo_ws = ax3.contourf(data['x'], data['y'], data['wspd_awo_qc'][0] - data['wspd_awo_ws_qc'][0], 
                        levels=wspd_diff_levels, cmap=diff_cmap, extend='both')


hcb = fig.colorbar(wspd_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax3, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=labelsize, pad=0.002)

add_storm_relative_essentials(ax3)
add_axis_labels(ax3, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax3, x_pos, y_pos, '(c)', fontsize=fontsize)
ax3.set_title('$AWO$-$CTL$ - $AWO_{ws}$-$EXP$ $U_{10}}$ ($m/s$)', fontsize=fontsize, pad=1)
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
add_axis_labels(ax4, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax4, x_pos, y_pos, '(d)', fontsize=fontsize)

ax4.set_title('$AWO$-$CTL$ $H_{s}}$ ($m$)', fontsize=fontsize, pad=1)
ax4.set_aspect('equal')

#AWO_ws
ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)

awo_ws_swh = ax5.contourf(data['x'], data['y'], data['swh_awo_ws_qc'][0], 
                        levels=swh_levels, cmap=swh_cmap, extend='both')

hcb = fig.colorbar(awo_ws_swh, shrink=shrink, aspect=aspect, ax=ax5, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax5)
add_axis_labels(ax5, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax5, x_pos, y_pos, '(e)', fontsize=fontsize)
ax5.set_title('$AWO_{ws}$-$EXP$ $H_{s}$ ($m$)', fontsize=fontsize, pad=1)
ax5.set_aspect('equal')



#Diff
ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)

swh_awo_awo_ws = ax6.contourf(data['x'], data['y'], data['swh_awo_qc'][0] - data['swh_awo_ws_qc'][0], 
                        levels=swh_diff_levels, cmap=diff_cmap, extend='both')

hcb = fig.colorbar(swh_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax6, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax6)
add_axis_labels(ax6, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax6, x_pos, y_pos, '(f)', fontsize=fontsize)
ax6.set_title('$AWO$-$CTL$ - $AWO_{ws}$-$EXP$ $H_{s}$ ($m$)', fontsize=fontsize, pad=1)
ax6.set_aspect('equal')

#Surface Temperatures
#AWO
ax7 = plt.subplot2grid(gridsize, (2, 0), colspan=1, rowspan=1)

awo_temp = ax7.contourf(data['x'], data['y'], data['sst_awo_qc'][0], 
                        cmap=sst_cmap,levels=sst_levels, extend='both')

hcb = fig.colorbar(awo_temp, shrink=shrink, aspect=aspect, ax=ax7, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax7)
add_axis_labels(ax7, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax7, x_pos, y_pos, '(g)', fontsize=fontsize)
ax7.set_title('$AWO$-$CTL$ $T$ ($^{\circ}C$)', fontsize=fontsize, pad=1)
ax7.set_aspect('equal')

#AWO_ws
ax8 = plt.subplot2grid(gridsize, (2, 1), colspan=1, rowspan=1)

awo_ws_temp = ax8.contourf(data['x'], data['y'], data['sst_awo_ws_qc'][0], 
                        cmap=sst_cmap,levels=sst_levels, extend='both')


hcb = fig.colorbar(awo_ws_temp, shrink=shrink, aspect=aspect, ax=ax8, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax8)
add_axis_labels(ax8, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax8, x_pos, y_pos, '(h)', fontsize=fontsize)
ax8.set_title('$AWO_{ws}$-$EXP$ $T$ ($^{\circ}C$)', fontsize=fontsize, pad=1)
ax8.set_aspect('equal')

#Diff
ax9 = plt.subplot2grid(gridsize, (2, 2), colspan=1, rowspan=1)

temp_awo_awo_ws = ax9.contourf(data['x'], data['y'], data['sst_awo_qc'][0] - data['sst_awo_ws_qc'][0], 
                        cmap=diff_cmap, levels=sst_diff_levels, extend='both')


hcb = fig.colorbar(temp_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax9, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax9)
add_axis_labels(ax9, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax9, x_pos, y_pos, '(i)', fontsize=fontsize)
ax9.set_title('$AWO$-$CTL$ - $AWO_{ws}$-$EXP$ $T$ ($^{\circ} C$)', fontsize=fontsize, pad=1)
ax9.set_aspect('equal')

#Surface Currents
#AWO
ax10 = plt.subplot2grid(gridsize, (3, 0), colspan=1, rowspan=1)

awo_sfc = ax10.contourf(data['x'], data['y'], data['sfc_awo_qc'][0], 
                         levels=sfc_levels, cmap=sfc_cmap, extend='max', zorder=zorder)


hcb = fig.colorbar(awo_sfc, shrink=shrink, aspect=aspect, ax=ax10, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax10)
add_axis_labels(ax10, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax10, x_pos, y_pos, '(j)', fontsize=fontsize)
ax10.set_title('$AWO$-$CTL$ $U$ ($m/s$)', fontsize=fontsize, pad=1)
ax10.set_aspect('equal')


#AWO_ws
ax11 = plt.subplot2grid(gridsize, (3, 1), colspan=1, rowspan=1)

awo_ws_sfc = ax11.contourf(data['x'], data['y'], data['sfc_awo_ws_qc'][0], 
                        levels=sfc_levels, cmap=sfc_cmap, extend='max')

hcb = fig.colorbar(awo_ws_sfc, shrink=shrink, aspect=aspect, ax=ax11, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax11)
add_axis_labels(ax11, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax11, x_pos, y_pos, '(k)', fontsize=fontsize)
ax11.set_title('$AWO_{ws}$-$EXP$ $U}$ ($m/s$)', fontsize=fontsize, pad=1)
ax11.set_aspect('equal')

# Diff
ax12 = plt.subplot2grid(gridsize, (3, 2), colspan=1, rowspan=1)

sfc_awo_awo_ws = ax12.contourf(data['x'], data['y'], data['sfc_awo_qc'][0] - data['sfc_awo_ws_qc'][0], 
                         levels=sfc_diff_levels, cmap=diff_cmap, extend='both', zorder=zorder)

hcb = fig.colorbar(sfc_awo_awo_ws, shrink=shrink, aspect=aspect, ax=ax12, pad=0.02)
hcb.ax.tick_params(color='k', length=1, width=0.5, labelsize=4, pad=0.002)

add_storm_relative_essentials(ax12)
add_axis_labels(ax12, fontsize, labelpad, labelsize, ticks[::2], ticks[::2])
add_corner_label(ax12, x_pos, y_pos, '(l)', fontsize=fontsize)
ax12.set_title('$AWO$-$CTL$ - $AWO_{ws}$-$EXP$ $U$ ($m/s$)', fontsize=fontsize, pad=1)
ax12.set_aspect('equal')



fig.tight_layout(pad=0.025, w_pad=0.025, h_pad=0.0025)

plt.savefig(PNG + 'storm_relative_average_wspd_swh_sst_sfc.py.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)
