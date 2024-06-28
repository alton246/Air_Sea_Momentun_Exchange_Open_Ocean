import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime 
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import xarray as xr
import cmaps
import pylab

from datetime import datetime,timedelta
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rcParams
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from helpers import *

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt


#Reading data
# Modelled Data
PATH_Hur = '/home/disk/orca/adaley17/Research/Stress_Separation/Hurricane_Earl/Data/'
PNG='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'
PNG2='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentum_Exchange_TC_Coast/Notebooks/Figure1/'

#UWIN-CM Best Track
earl_awo = PATH_Hur + 'earl_5.log_awo.csv'
earl_awo_ws = PATH_Hur + 'earl_5.log_awo-ws.csv'

#Ocean Data
ocn_file = 'archv.2010_244_01.nc'
awo_hyc_data = 'awo5_2010082700_gfs_3.7.1/' + ocn_file
awo_ws_hyc_data = 'awo5-ws_2010082700_gfs_3.7.1/' + ocn_file

#Accessing Ocean Data
awo_ocn_data  = xr.open_dataset(PATH_Hur + awo_hyc_data)
awo_ws_ocn_data  = xr.open_dataset(PATH_Hur + awo_ws_hyc_data)

#Computing Currents
awo_ocn_data['currents'] = np.sqrt(awo_ocn_data['u-vel']**2 + awo_ocn_data['v-vel']**2)
awo_ws_ocn_data['currents'] = np.sqrt(awo_ws_ocn_data['u-vel']**2 + awo_ws_ocn_data['v-vel']**2)

#Wave Data
wave_file = 'umwmout_2010-09-01_01:00:00.nc'
awo_umwm_data = 'awo5_2010082700_gfs_3.7.1/output/' + wave_file
awo_ws_umwm_data = 'awo5-ws_2010082700_gfs_3.7.1/output/' + wave_file

#Accessing Ocean Data
awo_wave_data  = xr.open_dataset(PATH_Hur + awo_umwm_data)
awo_ws_wave_data  = xr.open_dataset(PATH_Hur + awo_ws_umwm_data)

#Computing Ocean Surface Currents
awo_wave_data['currents'] = np.sqrt(awo_wave_data['uc']**2 + awo_wave_data['vc']**2)
awo_ws_wave_data['currents'] = np.sqrt(awo_ws_wave_data['uc']**2 + awo_ws_wave_data['vc']**2)


# IBTRAC Data
PATH_BT = '/home/orca/data/best_track/IBTrACS/'
track_file = 'ibtracs.NA.list.v04r00.csv'

#Reading track files
track_awo = atcf_csv(earl_awo) #AWO
track_awo_ws = atcf_csv(earl_awo_ws) #AWO-ws

#Dates of interest
date_of_int = pd.date_range('2010-08-27T00:00:00', '2010-09-03T00:00:00',
                                      freq='1H')

# UWINCM Storm Centers
#AWO
storm_centers_awo = extract_Storm_centers_over_period_time(date_of_int, track_awo)

#AWO_ws
storm_centers_awo_ws = extract_Storm_centers_over_period_time(date_of_int, track_awo_ws)

#Enter datetime string for the start and end time for UWINCM tracks
awo_start_date_index = track_awo.time[::6].flatten().tolist().index(datetime(2010,8,27,0))
awo_end_date_index = track_awo.time[::6].flatten().tolist().index(datetime(2010,9,2,0))

awo_swh_start_date_index = track_awo.time[::3].flatten().tolist().index(datetime(2010,8,30,6))
awo_swh_end_date_index = track_awo.time[::3].flatten().tolist().index(datetime(2010,9,1,3))

#Start of Plotting

#Set Domain Size and Cartopy Options here
min_lat = 14
max_lat = 30
min_lon = -80
max_lon = -60

crs = ccrs.PlateCarree()

plot_area = [max_lon + 360.0, min_lon + 360.0, min_lat, max_lat]

#Set Figure Options Here
skip_vec = 30
sst_levels = np.arange(28, 29.7, 0.1)
wspd_levels = np.arange(5,57.5,2.5)
diff_levels = np.arange(-0.3, 0.35, 0.05)
swh_levels = np.arange(0,17,1)
fontsize=8
headwidth=4
rcParams['xtick.major.pad']='0'
xlocator =4
ylocator=4
fontsize=10
markersize=2
lw=1
shrink=0.95

aspect=18
width=0.5
length=4

skip=24
labelsize=8
y_pos=0.92

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

gridsize = (3, 1)
fig = plt.figure(figsize=(8, 8))

#AWO Wave
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, projection=crs)

awo_swh = ax1.contourf(awo_wave_data['lon'][0], awo_wave_data['lat'][0], awo_wave_data['swh'][0], 
                        cmap=cmaps.NMCRef, levels=swh_levels, extend='max', transform=ccrs.PlateCarree())

#Track Info
ax1.scatter(storm_centers_awo[:,0][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
         linestyle='-', facecolors='white', edgecolors='k', marker='o', s=5)

#Wave Vector
vec_wave = ax1.quiver(awo_wave_data['lon'][0][::skip_vec, ::skip_vec], awo_wave_data['lat'][0][::skip_vec, ::skip_vec], 
                     u_wvd[::skip_vec, ::skip_vec], v_wvd[::skip_vec, ::skip_vec], 
                     angles='xy', scale_units='xy', scale=1, color='magenta', 
                     headwidth=headwidth, linewidths=0.5, edgecolors='magenta')

#StormDirection vector
vec_wave = ax1.quiver(storm_centers_awo[:,0][::3][awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_end_date_index], 
                     x_comp_storm_dir, y_comp_storm_dir, 
                     angles='xy', scale_units='xy', scale=0.4, color='white', 
                     headwidth=headwidth, linewidths=0.25, edgecolors='white')

#Wind Vector
vec_wind = ax1.quiver(awo_wave_data['lon'][0][::skip_vec, ::skip_vec], awo_wave_data['lat'][0][::skip_vec, ::skip_vec], 
                     u_wind[::skip_vec, ::skip_vec], v_wind[::skip_vec, ::skip_vec], 
                     angles='xy', scale_units='xy', scale=1, color='k', 
                     headwidth=headwidth, linewidths=0.5, edgecolors='k')

#Colorbar
hcb = fig.colorbar(awo_swh, shrink=shrink, aspect=aspect, ax=ax1, pad=0.02)
hcb.ax.tick_params(color='k', length=3, width=1.5, labelsize=labelsize, pad=0.002)



Cartopy_Features(ax1, fontsize, plot_area, 2, 2, 'k')
ax1.set_title('$AWO$-$CTL$ $H_{s}$ $(m)$', fontsize=fontsize, pad=1)
add_corner_label(ax1, y_pos, '(a)')

#AWO WSPD
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, projection=crs)
awo_wspd = ax2.contourf(awo_wave_data['lon'][0], awo_wave_data['lat'][0], awo_wave_data['wspd'][0], 
                        cmap=cmaps.MPL_gist_rainbow_r, levels=wspd_levels, extend='both', transform=ccrs.PlateCarree())

hcb = fig.colorbar(awo_wspd, shrink=shrink, aspect=aspect, ax=ax2, pad=0.02)
hcb.ax.tick_params(color='k', length=length, width=width, labelsize=labelsize, pad=0.002)
hcb.ax.minorticks_on()

Cartopy_Features(ax2, fontsize, plot_area, 2, 2, 'k')

#StormDirection vector
vec_wave = ax2.quiver(storm_centers_awo[:,0][::3][awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_end_date_index], 
                     x_comp_storm_dir, y_comp_storm_dir, 
                     angles='xy', scale_units='xy', scale=0.4, color='white', 
                     headwidth=headwidth, linewidths=0.25, edgecolors='white')

ax2.scatter(storm_centers_awo[:,0][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
         linestyle='-', facecolors='white', edgecolors='k', marker='o', s=5)

ax2.set_title('$AWO$-$CTL$ $U_{10}$ ($m/s$)', fontsize=fontsize, pad=1)

add_corner_label(ax2, y_pos, '(b)')

# AWO SST
ax3 = plt.subplot2grid(gridsize, (2, 0), colspan=1, rowspan=1, projection=crs)

awo_sst = ax3.contourf(awo_ocn_data['longitude'], awo_ocn_data['latitude'], awo_ocn_data['temp'][0][0], 
                        cmap=cmaps.MPL_jet, levels=sst_levels, extend='both', transform=ccrs.PlateCarree())

hcb = fig.colorbar(awo_sst, shrink=shrink, aspect=aspect, ax=ax3, pad=0.02)
hcb.ax.tick_params(color='k', length=length, width=width, labelsize=labelsize, pad=0.002)
hcb.ax.minorticks_on()

#StormDirection vector
vec_wave = ax3.quiver(storm_centers_awo[:,0][::3][awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_end_date_index], 
                     x_comp_storm_dir, y_comp_storm_dir, 
                     angles='xy', scale_units='xy', scale=0.4, color='white', 
                     headwidth=headwidth, linewidths=0.25, edgecolors='white')

ax3.scatter(storm_centers_awo[:,0][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
            storm_centers_awo[:,1][::3][awo_swh_start_date_index:awo_swh_end_date_index], 
         linestyle='-', facecolors='white', edgecolors='k', marker='o', s=5)

Cartopy_Features(ax3, fontsize, plot_area, 2, 2, 'k')
ax3.set_title('$AWO$-$CTL$ SST ($^{\circ}$C)', fontsize=fontsize, pad=1)

add_corner_label(ax3, y_pos, '(c)')


fig.tight_layout(pad=0, w_pad=0.25, h_pad=0)

plt.savefig(PNG + 'large_scale_swh_wspd_sst.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)

plt.savefig(PNG2 + 'large_scale_swh_wspd_sst.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)