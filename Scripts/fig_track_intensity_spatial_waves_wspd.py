import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime 
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import xarray as xr
import cmaps

from datetime import datetime,timedelta
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rcParams
from helpers import *

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

###Setting up PATHS
PATH_BT = '/home/orca/data/best_track/IBTrACS/'
PATH_Hur = '/home/disk/orca/adaley17/Research/Stress_Separation/Hurricane_Earl/Data/'
PNG='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'
PNG2='/home/disk/orca/adaley17/public_html/tmp/Presentations/AGU2024/Figures/'

earl_awo = PATH_Hur + 'earl_5.log_awo.csv'
earl_awo_ws = PATH_Hur + 'earl_5.log_awo-ws.csv'

#Ocean Data
ocn_file = 'archv.2010_244_01.nc'
awo_hyc_data = 'awo5_2010082700_gfs_3.7.1/' + ocn_file
awo_ws_hyc_data = 'awo5-ws_2010082700_gfs_3.7.1/' + ocn_file

#Loading Ocean Data
awo_ocn_data  = xr.open_dataset(PATH_Hur + awo_hyc_data)
awo_ws_ocn_data  = xr.open_dataset(PATH_Hur + awo_ws_hyc_data)

# IBTRAC File
track_file = 'ibtracs.NA.list.v04r01.csv'

#Wave Data
wave_file = 'umwmout_2010-09-01_00:00:00.nc'
awo_umwm_data = 'awo5_2010082700_gfs_3.7.1/output/' + wave_file
# awo_ws_umwm_data = 'awo5-ws_2010082700_gfs_3.7.1/output/' + wave_file
awo_wave_data  = xr.open_dataset(PATH_Hur + awo_umwm_data)

#Accessing Wave Data
awo_wave_data  = xr.open_dataset(PATH_Hur + awo_umwm_data)
# awo_ws_wave_data  = xr.open_dataset(PATH_Hur + awo_ws_umwm_data)

#Reading trak files
track_awo = atcf_csv(earl_awo) #AWO
track_awo_ws = atcf_csv(earl_awo_ws) #AWO-ws

#Dates of interest
date_of_int = pd.date_range('2010-08-27T00:00:00', '2010-09-02T23:00:00',
                                      freq='1H')

# UWINCM Storm Centers
#AWO
storm_centers_awo = np.zeros((len(date_of_int),2), dtype=float)
for i in range(len(date_of_int)):
    lats_lons_awo = getStormCenter(date_of_int[i],track_awo)
    storm_centers_awo[i,0] = lats_lons_awo[0]
    storm_centers_awo[i,1] = lats_lons_awo[1]

#AWoO_ws
storm_centers_awo_ws = np.zeros((len(date_of_int),2), dtype=float)
for i in range(len(date_of_int)):
    lats_lons_awo_ws = getStormCenter(date_of_int[i],track_awo_ws)
    storm_centers_awo_ws[i,0] = lats_lons_awo_ws[0]
    storm_centers_awo_ws[i,1] = lats_lons_awo_ws[1]


#Enter storm name and year here
storm_name = 'EARL'
year_of_storm = '2010'

# Extrating the number of entries and data for Hurricane EARL
num_entries_earl, bt_data = counter(PATH_BT, track_file, storm_name, year_of_storm)

# Extracting pressure, wspd, date, lat, lon and 
# colors for Hurricane EARL from best track dataset

earl_bt_data = extract_ibt_data(num_entries_earl, bt_data, storm_name, year_of_storm)

#Enter datetime string for the start and end time for NHC best track
bt_start_date_index = earl_bt_data.date[::2].flatten().tolist().index(datetime(2010,8,27,0))
bt_end_date_index = earl_bt_data.date[::2].flatten().tolist().index(datetime(2010,9,2,6))

#Enter datetime string for the start and end time for UWINCM tracks
awo_start_date_index = track_awo.time[::6].flatten().tolist().index(datetime(2010,8,27,0))
awo_end_date_index = track_awo.time[::6].flatten().tolist().index(datetime(2010,9,2,0))

awo_swh_start_date_index = track_awo.time[::3].flatten().tolist().index(datetime(2010,8,30,6))
awo_swh_end_date_index = track_awo.time[::3].flatten().tolist().index(datetime(2010,9,1,3))

##Start of Plotting
#Cartopy dependencies
crs = ccrs.PlateCarree()
coast = cfeature.GSHHSFeature(scale='high', levels=[1,], edgecolor='k')
lakes = cfeature.GSHHSFeature(scale='high', levels=[2,], edgecolor='face',facecolor='grey')

#Track Domain
min_lat_track = 14
max_lat_track = 30
min_lon_track = -80
max_lon_track = -38
plot_area_track = [max_lon_track + 360.0, min_lon_track + 360.0, min_lat_track,max_lat_track]

#Spatial Plor Domain
min_lat = 14
max_lat = 30
min_lon = -80
max_lon = -60
plot_area = [max_lon + 360.0, min_lon + 360.0, min_lat, max_lat]

# #Figure Options
fontlabel_size=8
fonttick_size =8
wind_yticks=np.arange(10,80,10)
mslp_yticks = np.arange(930,1020, 10)
markersize=3
skip=24
fontsize=6
lw=0.75
size=12
bt_skip=2
mod_skip=6

#Colorbar Options
sst_levels = np.arange(28, 29.7, 0.1) #SSt levels
swh_levels = np.arange(0,17,1) #SWH levels
rcParams['xtick.major.pad']='0'
shrink=0.95
aspect=18
width=0.5
length=4
labelsize=8

#Vector Options
skip_vec = 30
headwidth=4

#Figure Label Options
x_pos=0.05
y_pos=0.86

#Wave Vector Components
u_wvd = np.cos(awo_wave_data['mwd'][0])
v_wvd = np.sin(awo_wave_data['mwd'][0])

#Wind Vector components
u_wind = np.cos(awo_wave_data['wdir'][0])
v_wind = np.sin(awo_wave_data['wdir'][0])

#Storm Direction Vector Components
awo_storm_dir = getStormDirection(plain2datetime('2010090101'), track_awo)
x_comp_storm_dir = np.cos(awo_storm_dir)
y_comp_storm_dir = np.sin(awo_storm_dir)

# print(awo_ocn_data['temp'][0][0].values)

gridsize = (2, 4)
fig = plt.figure(figsize=(5,3))

#Track
ax1 = plt.subplot2grid(gridsize, (0,0), colspan=2, rowspan=1, projection=crs)


Cartopy_Features(ax1, 8, plot_area_track, 4, 4, 'k')
#Scatter
ax1.scatter(earl_bt_data.lons[::bt_skip], earl_bt_data.lats[::bt_skip], marker="s", color='k',
            edgecolors='k', s=size, label='best track')
ax1.scatter(storm_centers_awo[:,0][::mod_skip], storm_centers_awo[:,1][::mod_skip], marker="o", 
            facecolors='red',s=size)
ax1.scatter(storm_centers_awo_ws[:,0][::mod_skip], storm_centers_awo_ws[:,1][::mod_skip], marker="*", 
            facecolors='cyan', s=size)

#Line Plots
ax1.plot(earl_bt_data.lons, earl_bt_data.lats, linestyle='-', color='k', linewidth=lw, label='best track')
ax1.plot(storm_centers_awo[:,0], storm_centers_awo[:,1], linewidth=lw, linestyle='-', color='red', label='$CTL$')
ax1.plot(storm_centers_awo_ws[:,0], storm_centers_awo_ws[:,1], linewidth=lw, linestyle='-', color='cyan', label='$EXP$')


for i in range(len(track_awo.time[::skip])):
    print(track_awo.time[::skip][i])

    ax1.text(storm_centers_awo[:,0][::skip][i]+1, storm_centers_awo[:,1][::skip][i]+0.8, 
        track_awo.time[::skip][i].strftime('%m-%d'), fontsize=5.5, fontweight='semibold')
    x_awo = [storm_centers_awo[:,0][::skip][i], storm_centers_awo[:,0][::skip][i]+1]
    y_awo = [storm_centers_awo[:,1][::skip][i], storm_centers_awo[:,1][::skip][i]+1]

    x_awo_ws = [storm_centers_awo_ws[:,0][::skip][i], storm_centers_awo[:,0][::skip][i]+1]
    y_awo_ws = [storm_centers_awo_ws[:,1][::skip][i], storm_centers_awo[:,1][::skip][i]+1]
    # print(storm_centers_awo[:,0][::skip])
    x_bt = [earl_bt_data.lons[::8][3:10][i][0], storm_centers_awo[:,0][::skip][i]+1]
    y_bt = [earl_bt_data.lats[::8][3:10][i][0], storm_centers_awo[:,1][::skip][i]+1]
    ax1.plot(x_awo, y_awo, 'red', linestyle="-", lw=lw)
    ax1.plot(x_awo_ws, y_awo_ws, 'cyan', linestyle="-", lw=lw)
    # print(x_bt)
    ax1.plot(x_bt, y_bt, 'k', linestyle="-", lw=0.75) #Needs to be corrected

#Legend
# Track_Legend(axis, fontsize, markersize, xlocator, ylocator, loc)
Track_Legend_v2(ax1, 5, 5, 0.84, 0.78, 'center')
# plt.legend(loc='upper left', prop = {"size": 5}, frameon=True, fancybox=True, shadow=True)
add_corner_label(ax1, x_pos, y_pos, '(a)', fontsize=5.5)


# Intensisty
# WSPD
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=2, rowspan=1)
ax2.plot(earl_bt_data.date[::bt_skip][bt_start_date_index:bt_end_date_index], earl_bt_data.max_ws[::bt_skip][bt_start_date_index:bt_end_date_index], 
         color='k', markersize=markersize, linewidth=lw, linestyle='-', marker='s', label='best track')
ax2.plot(track_awo.time[::mod_skip], track_awo.wspd[::mod_skip], linewidth=lw, color='red', linestyle='-', 
         marker='o', markersize=markersize, label='$CTL$')
ax2.plot(track_awo_ws.time[::mod_skip], track_awo_ws.wspd[::mod_skip], color='cyan', linestyle='-', 
         marker='*', linewidth=lw, markersize=markersize, label='$EXP$')
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d")) #Date Format
ax2.set_xlim([track_awo.time[::mod_skip][awo_start_date_index], 
              track_awo.time[::mod_skip][awo_end_date_index]]) #Minimum and Max x values 
# ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d")) #Date Format
ax2.set_ylabel('Wind Speed (m/s)', fontsize=fontlabel_size)
ax2.set_xlabel('Date (mm-dd)', fontsize=fontlabel_size)
add_corner_label(ax2, x_pos, y_pos, '(c)', fontsize=5.5)
ax2.set_yticks(ticks=wind_yticks)
ax2.set_xticks(ticks = track_awo.time[::24])
ax2.set_ylim([15, 70] ) #Minimum and Max x values 
ax2.set_xlim([track_awo.time[::mod_skip][awo_start_date_index], 
              track_awo.time[::mod_skip][awo_end_date_index]]) #Minimum and Max x values 
ax2.tick_params(axis='both', which='major', labelsize=fontlabel_size)
ax2.legend(loc='lower right', prop = {"size": 5}, frameon=True, fancybox=True, shadow=True) #Legend
ax2.grid(linestyle=':')

#Waves
#AWO Wave
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=2, rowspan=1, projection=crs)

awo_swh = ax3.contourf(awo_wave_data['lon'][0], awo_wave_data['lat'][0], awo_wave_data['swh'][0], 
                        cmap=cmaps.NMCRef, levels=swh_levels, extend='max', transform=ccrs.PlateCarree())

#Track Info
ax3.scatter(storm_centers_awo[:,0][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
         linestyle='-', facecolors='white', edgecolors='k', marker='o', s=5)

#Wave Vector
vec_wave = ax3.quiver(awo_wave_data['lon'][0][::skip_vec, ::skip_vec], awo_wave_data['lat'][0][::skip_vec, ::skip_vec], 
                     u_wvd[::skip_vec, ::skip_vec], v_wvd[::skip_vec, ::skip_vec], 
                     angles='xy', scale_units='xy', scale=1, color='magenta', 
                     headwidth=headwidth, linewidths=0.5, edgecolors='magenta')

#StormDirection vector
vec_wave = ax3.quiver(storm_centers_awo[:,0][::3][awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_end_date_index], 
                     x_comp_storm_dir, y_comp_storm_dir, 
                     angles='xy', scale_units='xy', scale=0.4, color='white', 
                     headwidth=headwidth, linewidths=0.25, edgecolors='white')

#Wind Vector
vec_wind = ax3.quiver(awo_wave_data['lon'][0][::skip_vec, ::skip_vec], awo_wave_data['lat'][0][::skip_vec, ::skip_vec], 
                     u_wind[::skip_vec, ::skip_vec], v_wind[::skip_vec, ::skip_vec], 
                     angles='xy', scale_units='xy', scale=1, color='k', 
                     headwidth=headwidth, linewidths=0.5, edgecolors='k')

#Colorbar
hcb = fig.colorbar(awo_swh, shrink=shrink, aspect=aspect, ax=ax3, pad=0.02)
hcb.ax.tick_params(color='k', length=3, width=1.5, labelsize=labelsize, pad=0.002)



Cartopy_Features(ax3, fontsize, plot_area, 2, 2, 'k')
ax3.set_title('$CTL$ $H_{s}$ $(m)$', fontsize=fontsize, pad=1)
add_corner_label(ax3, x_pos, y_pos, '(b)', 5.5)

ax4 = plt.subplot2grid(gridsize, (1, 2), colspan=2, rowspan=1, projection=crs)

awo_sst = ax4.contourf(awo_ocn_data['longitude'], awo_ocn_data['latitude'], awo_ocn_data['temp'][0][0], 
                        cmap=cmaps.MPL_jet, levels=sst_levels, extend='both', transform=ccrs.PlateCarree())

hcb = fig.colorbar(awo_sst, shrink=shrink, aspect=aspect, ax=ax4, pad=0.02)
hcb.ax.tick_params(color='k', length=length, width=width, labelsize=labelsize, pad=0.002)
hcb.ax.minorticks_on()

#StormDirection vector
vec_wave = ax4.quiver(storm_centers_awo[:,0][::3][awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_end_date_index], 
                     x_comp_storm_dir, y_comp_storm_dir, 
                     angles='xy', scale_units='xy', scale=0.4, color='white', 
                     headwidth=headwidth, linewidths=0.25, edgecolors='white')

ax4.scatter(storm_centers_awo[:,0][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
         linestyle='-', facecolors='white', edgecolors='k', marker='o', s=5)

Cartopy_Features(ax4, fontsize, plot_area, 2, 2, 'k')
ax4.set_title('$CTL$ SST ($^{\circ}$C)', fontsize=fontsize, pad=1)

add_corner_label(ax4, x_pos, y_pos, '(d)', fontsize)


fig.tight_layout(pad=0, w_pad=0.5, h_pad=0)


plt.savefig(PNG + 'track_wspd_intensity_spatial_waves_sst.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)