import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime 
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import xarray as xr
import cmaps
import os
import glob

from datetime import datetime,timedelta
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rcParams
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from helpers import *
from scipy.spatial import cKDTree

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
from concurrent.futures import ThreadPoolExecutor
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt


# Modelled Data
PATH_Hur = '/home/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/'
DWM_PATH = '/home/disk/oldhot/from_orca/from_miami/uwin2/milan/output/earl/'
NDBC_PATH = '/home/disk/orca/adaley17/orca3/adaley17/data/ndbc/'
PNG='/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'

#Buoy Files 
buoy_file_43 = '41043h2010.txt' #41043
buoy_file_46 = '41046h2010.txt' #41046

#Accessing Buoy Data
df_043 = pd.read_csv(NDBC_PATH + buoy_file_43, delim_whitespace=True, skiprows=[1], header=0)
df_046 = pd.read_csv(NDBC_PATH + buoy_file_46, delim_whitespace=True, skiprows=[1], header=0)

#UWIN-CM Best Track
earl_awo = PATH_Hur + 'earl_5.log_awo.csv'
earl_awo_ws = PATH_Hur + 'earl_5.log_awo-ws.csv'

#Reading track files
track_awo = atcf_csv(earl_awo) #AWO
track_awo_ws = atcf_csv(earl_awo_ws) #AWO-ws

#Dates of interest
date_of_int = pd.date_range('2010-08-27T00:00:00', '2010-09-02T23:00:00',
                                      freq='1h')

# UWINCM Storm Centers
#AWO
storm_centers_awo = np.zeros((len(date_of_int),2), dtype=float)
for i in range(len(date_of_int)):
    lats_lons_awo = getStormCenter(date_of_int[i],track_awo)
    storm_centers_awo[i,0] = lats_lons_awo[0]
    storm_centers_awo[i,1] = lats_lons_awo[1]

#Buoy Locations
buoy_names = ['41043', '41046', '41044']
buoy_lat = [21.026, 23.822, 21.582]
buoy_lon = [-64.793, -68.393, -58.63]

#Comparing UWINCM SWH to Buoy Data

os.chdir(DWM_PATH + 'awo5_2010082700_gfs_3.7.1/') #Getting into the right directory
wave_files = sorted(glob.glob("umwmout*.nc")) #Grabbing all the wave files in the directory

#Extracting UWIN-CM SWH time series for each buoy
awo_swh_046, awo_ws_swh_046 = extract_swh_series(wave_files, DWM_PATH, 'swh', buoy_lat[1], buoy_lon[1])
awo_swh_043, awo_ws_swh_043 = extract_swh_series(wave_files, DWM_PATH, 'swh', buoy_lat[0], buoy_lon[0])


buoy_043_qc = clean_buoy_data(df_043) #Cleaning Buoy Data
buoy_046_qc = clean_buoy_data(df_046) #Cleaning Buoy Data

subset_df_043 = buoy_043_qc['2010-08-26 23:50:00':'2010-09-01 11:50:00'] #Subsetting Buoy Data
subset_df_046 = buoy_046_qc['2010-08-26 23:50:00':'2010-09-01 11:50:00']#Subsetting Buoy Data

buoy_date_int = pd.date_range('2010-08-27T00:00:00', '2010-09-01T12:00:00',
                                      freq='1h')
start_index = 0
end_index= 133

# Figure Options
xlocator =4
ylocator=4
fontsize=6
markersize=2
lw=1
shrink=0.50

aspect=18
width=0.5
length=4

skip_track=6
skip_date=24
labelsize=6
x_pos=0.04
y_pos=0.88

#Start of Plotting

#Set Domain Size and Cartopy Options here
min_lat = 14
max_lat = 30
min_lon = -80
max_lon = -54

ymin, ymax = 0, 15
yticks= np.arange(0,18,2)

crs = ccrs.PlateCarree()

plot_area = [max_lon + 360.0, min_lon + 360.0, min_lat, max_lat]

gridsize = (2, 2)
fig = plt.figure(figsize=(4, 4))

#AWO Wave
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=2, rowspan=1, projection=crs)

ax1.scatter(track_awo.lon[::skip_track], track_awo.lat[::skip_track], c='k', s=markersize)
ax1.scatter(track_awo_ws.lon[::skip_track], track_awo_ws.lat[::skip_track], c='r', s=markersize)
ax1.scatter(buoy_lon, buoy_lat, c='b', s=markersize*2, marker='*')

#Add Buoy names
for name, lon, lat in zip(buoy_names, buoy_lon, buoy_lat):
    ax1.text(lon, lat + 0.5, name, color='b', fontsize=fontsize, ha='center', va='bottom')


for i in range(-4, 0):
    print(track_awo.time[::skip_date][i])

    ax1.text(storm_centers_awo[:,0][::skip_date][i]+0.2, storm_centers_awo[:,1][::skip_date][i]+0.2, 
        track_awo.time[::skip_date][i].strftime('%m-%d'), fontsize=6, fontweight='semibold')


Cartopy_Features(ax1, fontsize, plot_area, 2, 2, 'k')
# ax1.set_title('$CTL$ $H_{s}$ $(m)$', fontsize=fontsize, pad=1)
add_corner_label(ax1, x_pos, y_pos, '(a)', fontsize)


ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)
ax2.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%m-%d'))
ax2.plot(subset_df_043.index, subset_df_043['WVHT'], color='k', label='Buoy 41043')
ax2.plot(buoy_date_int, awo_swh_043[start_index:end_index], color='red', label='CTL')
ax2.plot(buoy_date_int, awo_ws_swh_043[start_index:end_index], color='cyan', label='EXP')

ax2.set_ylim(ymin, ymax)


ax2.set_xlabel('Date (mm-dd)', fontsize=fontsize)
ax2.set_ylabel('Sig. Wave Height (m)', fontsize=fontsize)
ax2.tick_params(axis='both', labelsize=fontsize)
ax2.set_yticks(yticks)
ax2.grid(True, linestyle=':', linewidth=1)
ax2.legend(loc='center left', fontsize=fontsize, shadow=True)
add_corner_label(ax2, x_pos, y_pos, '(b)', fontsize)


ax3 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)
ax3.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%m-%d'))
ax3.plot(subset_df_046.index, subset_df_046['WVHT'], color='k', label='Buoy 41046', )
ax3.plot(buoy_date_int, awo_swh_046[start_index:end_index], color='red', label='CTL')
ax3.plot(buoy_date_int, awo_ws_swh_046[start_index:end_index], color='cyan', label='EXP')

ax3.set_ylim(ymin, ymax)
ax3.set_xlabel('Date mm-dd', fontsize=fontsize)
ax3.set_ylabel('Sig. Wave Height (m)', fontsize=fontsize)
ax3.tick_params(axis='both', labelsize=fontsize)
ax3.set_yticks(yticks)
ax3.grid(True, linestyle=':', linewidth=1)
ax3.legend(loc='center left', fontsize=fontsize, shadow=True)
add_corner_label(ax3, x_pos, y_pos, '(c)', fontsize)



fig.tight_layout(pad=0, w_pad=0.25, h_pad=1)

plt.savefig(PNG + 'buoy_wave_analysis.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)
