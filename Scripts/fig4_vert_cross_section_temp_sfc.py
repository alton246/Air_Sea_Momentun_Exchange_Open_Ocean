import xarray as xr 
import numpy as np
import pandas as pd
import cmaps
from helpers import *

#Setting PATH of data and 
PATH = '/home/disk/orca/adaley17/Research/Stress_Separation/Hurricane_Earl/Notebooks/SST_Analysis/Composites/Data_Currents_Zlevels/'
PNG='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'


#Opening data files
with xr.open_mfdataset(PATH + 'storm_relative_hycom_3d.20100901??.nc', use_cftime=True) as ds:
    ds_timemean_slice = ds.copy().sel(y=-98, method='nearest').sel(time=slice('2010-09-01T00:00:00','2010-09-01T05:00:00')).mean(dim='time')

# Accessing Data of interest
X = ds_timemean_slice['x'].data
Y = ds_timemean_slice['z_interp'].data
awo_temp = ds_timemean_slice['temp_interp_awo'].data
awo_ws_temp = ds_timemean_slice['temp_interp_awo_ws'].data

#Computing Currents
awo_curr = np.sqrt(ds_timemean_slice['u-vel_interp_awo']**2 + ds_timemean_slice['v-vel_interp_awo']**2)
awo_ws_curr = np.sqrt(ds_timemean_slice['u-vel_interp_awo_ws']**2 + ds_timemean_slice['v-vel_interp_awo_ws']**2)

#Contour levels
curr_levels = np.arange(0,2.10,0.10)
sst_levels=np.arange(25, 30.25, 0.25)
diff_levels=np.arange(-0.2,0.25,0.05)
sst_diff_levels = np.arange(-0.6, 0.7, 0.1)

#Colorbar Options
size=0.25
shrink=0.25
cur_cmaps = cmaps.BkBlAqGrYeOrReViWh200_r
sst_cmaps = cmaps.MPL_jet

# Set Font information
fontsize=5.5
labelsize=5.5
labelpad=0.10

#Set tick size here
yticks=np.arange(0,125,25)
xticks=np.arange(-100, 120, 20)

#Set x and y pos here
x_pos = 0.065
y_pos = 0.8

gridsize = (2, 3)
fig = plt.figure(figsize=(6, 4))

#Currents
#AWO
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)

temp_awo = ax1.contourf(X, Y, awo_temp,
                  levels=sst_levels, extend='both', cmap=sst_cmaps)

#Colorbar
hcb = fig.colorbar(temp_awo, shrink=shrink, aspect=10, ax=ax1, pad=0.02)
hcb.ax.tick_params(color='k', length=2.0, width=1.0, labelsize=5, pad=0.002)

#Mix laye depth
# ax1.scatter(awo_x_dist, awo_vcross_data['mix_dpth'][0]/9800,
#             c='w', s=size, edgecolors='k')

ax1.set_aspect('equal')
ax1.set_xlim([-100,100])
ax1.set_ylim([100,1])
cross_section_essential(30, 90, 10, 0.5)
add_vertical_cross_section_axis_labels(ax1, fontsize, labelpad, labelsize, xticks[::2], yticks)
ax1.set_title('$AWO$-$CTL$  $T$ ($^{\circ}C$)', fontsize=fontsize, pad=1)

add_corner_label(ax1, x_pos, y_pos,'(a)', fontsize)

#AWO_ws
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)

temp_awo_ws = ax2.contourf(X, Y, awo_ws_temp,
                  levels=sst_levels, extend='both', cmap=sst_cmaps)

hcb = fig.colorbar(temp_awo_ws,shrink=shrink, aspect=10, ax=ax2, pad=0.02)
hcb.ax.tick_params(color='k', length=2.0, width=1.0, labelsize=5, pad=0.002)


ax2.set_aspect('equal')
ax2.set_xlim([-100,100])
ax2.set_ylim([100,1])
cross_section_essential(30, 90, 10, 0.5)
add_vertical_cross_section_axis_labels(ax2, fontsize, labelpad, labelsize, xticks[::2], yticks)
ax2.set_title(' $AWO_{ws}$-$EXP$ $T$ ($^{\circ}C$)', fontsize=fontsize, pad=1)

add_corner_label(ax2, x_pos, y_pos,'(b)', fontsize)

ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)

# ax3.pcolormesh(awo_X, awo_Y, awo_vcross_data['curr_speed'][0] - awo_ws_vcross_data['curr_speed'][0],
#                   vmin=-0.5, vmax=0.5, cmap=cmaps.MPL_bwr, rasterized=True)
temp_diff_diff = ax3.contourf(X, Y, awo_temp - awo_ws_temp,
                  levels=sst_diff_levels, extend='both', cmap=cmaps.MPL_bwr)

# plt.clabel(temp_diff, inline = True, 
#            fontsize=fontsize, fmt='%1.2f', colors = 'k')

#Colorbar
hcb = fig.colorbar(temp_diff_diff, shrink=shrink, aspect=10, ax=ax3, pad=0.02)
hcb.ax.tick_params(color='k', length=2.0, width=1.0, labelsize=5, pad=0.002)


ax3.set_aspect('equal')
ax3.set_xlim([-100,100])
ax3.set_ylim([100,1]) # Limits
cross_section_essential(30, 75, 10, 0.5)
add_vertical_cross_section_axis_labels(ax3, fontsize, labelpad, labelsize, xticks[::2], yticks)
ax3.set_aspect('equal')

ax3.set_title('$AWO$-$CTL$ $ - $ $AWO_{ws}$-$EXP$ $T$ ($^{\circ}C$)', loc='center', fontsize=fontsize, pad=1)
add_corner_label(ax3, x_pos, y_pos,'(c)', fontsize)

ax4 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)

curr_awo = ax4.contourf(X, Y, awo_curr,
                  levels=curr_levels,extend='max', cmap=cur_cmaps)

#Colorbar
hcb = fig.colorbar(curr_awo,shrink=shrink, aspect=10, ax=ax4, pad=0.02)
hcb.ax.tick_params(color='k', length=2.0, width=1.0, labelsize=5, pad=0.002)

cross_section_essential(30, 75, 10, 0.5)
add_vertical_cross_section_axis_labels(ax4, fontsize, labelpad, labelsize, xticks[::2], yticks)
ax4.set_xlim([-100,100])
ax4.set_ylim([100,0]) # Limits
ax4.set_aspect('equal')

ax4.set_title('$AWO$-$CTL$ $U$ ($m/s$)', fontsize=fontsize, pad=1)
add_corner_label(ax4, x_pos, y_pos,'(d)', fontsize)



ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)

curr_awo_ws = ax5.contourf(X, Y, awo_curr,
                  levels=curr_levels,extend='max', cmap=cur_cmaps)

#Colorbar
hcb = fig.colorbar(curr_awo_ws,shrink=shrink, aspect=10, ax=ax5, pad=0.02)
hcb.ax.tick_params(color='k', length=2.0, width=1.0, labelsize=5, pad=0.002)

cross_section_essential(30, 75, 10, 0.5)
add_vertical_cross_section_axis_labels(ax5, fontsize, labelpad, labelsize, xticks[::2], yticks)

ax5.set_aspect('equal')
ax5.set_ylim([100,0]) # Limits
ax5.set_xlim([-100,100])
ax5.set_title(' $AWO_{ws}$-$EXP$ $U$ ($m/s$)', fontsize=fontsize, pad=1)
add_corner_label(ax5, x_pos, y_pos,'(e)', fontsize)


ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)

curr_diff = ax6.contourf(X, Y, awo_curr - awo_ws_curr,
                  levels=diff_levels,extend='both', cmap=cmaps.MPL_bwr)

#Colorbar
hcb = fig.colorbar(curr_diff,shrink=shrink, aspect=10, ax=ax6, pad=0.02)
hcb.ax.tick_params(color='k', length=2.0, width=1.0, labelsize=5, pad=0.002)

# #Mix Layer Depth
# ax6.scatter(awo_x_dist, awo_vcross_data['mix_dpth'][0]/9800,
#             c='cyan', s=size, edgecolors='cyan', label='$AWO$-$CTL$')
# ax6.scatter(awo_ws_x_dist, awo_ws_vcross_data['mix_dpth'][0]/9800,
#             c='k', s=size,edgecolors='k', label='$AWO_{ws}$-$EXP$')

cross_section_essential(30, 75, 10, 0.5)
add_vertical_cross_section_axis_labels(ax6, fontsize, labelpad, labelsize, xticks[::2], yticks)

ax6.set_aspect('equal')
ax6.set_xlim([-100,100])
ax6.set_ylim([100,0]) # Limits
ax6.set_title('$AWO$-$CTL$ $ - $ $AWO_{ws}$-$EXP$ \n$U$ ($m/s$)', loc='center', fontsize=fontsize, pad=1)
add_corner_label(ax6, x_pos, y_pos,'(f)', fontsize)


fig.tight_layout(pad=0.5, w_pad=0.5, h_pad=-8)

plt.savefig(PNG +'fig4_vert_cross_sfc_temp_2010090100-2010090105.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)