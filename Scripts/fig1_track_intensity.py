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
from helpers import *

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

###Setting up PATHS

PATH_Hur = '/home/disk/orca/adaley17/Research/Stress_Separation/Hurricane_Earl/Data/'
PNG='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'
PNG2='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentum_Exchange_TC_Coast/Notebooks/Figure1/'

earl_awo = PATH_Hur + 'earl_5.log_awo.csv'
earl_awo_ws = PATH_Hur + 'earl_5.log_awo-ws.csv'

# IBTRAC Data
PATH_BT = '/home/orca/data/best_track/IBTrACS/'
track_file = 'ibtracs.NA.list.v04r00.csv'

#Wave Data
wave_file = 'umwmout_2010-09-01_00:00:00.nc'
awo_umwm_data = 'awo5_2010082700_gfs_3.7.1/output/' + wave_file
# awo_ws_umwm_data = 'awo5-ws_2010082700_gfs_3.7.1/output/' + wave_file
awo_wave_data  = xr.open_dataset(PATH_Hur + awo_umwm_data)

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

earl_bt_data = extract_ibt_data(num_entries_earl, bt_data)

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

# #Figure Options
fontlabel_size=8
fonttick_size =8
wind_yticks=np.arange(10,80,10)
mslp_yticks = np.arange(930,1020, 10)
markersize=3
skip=24
fontsize=10
lw=0.75
size=12
bt_skip=2
mod_skip=6


gridsize = (2, 1)
fig = plt.figure(figsize=(4,4))

#Track
ax1 = plt.subplot2grid(gridsize, (0,0), colspan=1, rowspan=1, projection=crs)


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
ax1.plot(storm_centers_awo[:,0], storm_centers_awo[:,1], linewidth=lw, linestyle='-', color='red', label='$AWO$')
ax1.plot(storm_centers_awo_ws[:,0], storm_centers_awo_ws[:,1], linewidth=lw, linestyle='-', color='cyan', label='$AWO_{ws}$')

for i in range(len(track_awo.time[::skip])):
    print(track_awo.time[::skip][i])

    ax1.text(storm_centers_awo[:,0][::skip][i]+1, storm_centers_awo[:,1][::skip][i]+0.8, 
        track_awo.time[::skip][i].strftime('%m-%d'), fontsize=6, fontweight='semibold')
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
Track_Legend_v2(ax1, 5, 5, 0.88, 0.85, 'center')
# plt.legend(loc='upper left', prop = {"size": 5}, frameon=True, fancybox=True, shadow=True)
add_corner_label(ax1, '(a)')


# Intensisty
# WSPD
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)
ax2.plot(earl_bt_data.date[::bt_skip][bt_start_date_index:bt_end_date_index], earl_bt_data.max_ws[::bt_skip][bt_start_date_index:bt_end_date_index], 
         color='k', markersize=markersize, linewidth=lw, linestyle='-', marker='s', label='best track')
ax2.plot(track_awo.time[::mod_skip], track_awo.wspd[::mod_skip], linewidth=lw, color='red', linestyle='-', 
         marker='o', markersize=markersize, label='$AWO$-$CTL$')
ax2.plot(track_awo_ws.time[::mod_skip], track_awo_ws.wspd[::mod_skip], color='cyan', linestyle='-', 
         marker='*', linewidth=lw, markersize=markersize, label='$AWO_{ws}$-$EXP$')
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d")) #Date Format
ax2.set_xlim([track_awo.time[::mod_skip][awo_start_date_index], 
              track_awo.time[::mod_skip][awo_end_date_index]]) #Minimum and Max x values 
# ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d")) #Date Format
ax2.set_ylabel('Wind Speed (m/s)', fontsize=fontlabel_size)
ax2.set_xlabel('Date (mm-dd)', fontsize=fontlabel_size)
add_corner_label(ax2, '(b)')
ax2.set_yticks(ticks=wind_yticks)
ax2.set_xticks(ticks = track_awo.time[::24])
ax2.set_ylim([15, 70] ) #Minimum and Max x values 
ax2.set_xlim([track_awo.time[::mod_skip][awo_start_date_index], 
              track_awo.time[::mod_skip][awo_end_date_index]]) #Minimum and Max x values 
ax2.tick_params(axis='both', which='major', labelsize=fontlabel_size)
ax2.legend(loc='lower right', prop = {"size": 5}, frameon=True, fancybox=True, shadow=True) #Legend
ax2.grid(linestyle=':')


plt.savefig(PNG2 + 'track_wspd_intensity.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)