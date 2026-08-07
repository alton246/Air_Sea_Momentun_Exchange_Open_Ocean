import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import glob
import cmaps
import warnings
warnings.filterwarnings('ignore')

## Import Cartopy stuff.
import cartopy.crs as ccrs
from scipy.spatial import cKDTree
from helpers import *
from math import radians, sin, cos, sqrt, atan2
from scipy.ndimage import gaussian_filter


# PNG='/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'

#Model data directory
Model_DIR = '/home/orca3/adaley17/Projects/air_sea_mom_exchange_open_ocean/'
PNG='/home/disk/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'

#Obs data directory
Obs_DIR = '/home/disk/orca3/adaley17/data/sfmr/'

#Reconaissance SFMR data file
sfmr_file = 'AFRC_SFMR20100901U1.nc'

#Track data
awo_track_file = 'earl_5.log_awo.csv'
awo_ws_track_file = 'earl_5.log_awo-ws.csv'

#Load the track data
awo_track = atcf_csv(Model_DIR + awo_track_file)
awo_ws_track = atcf_csv(Model_DIR + awo_ws_track_file)

#Load the SFMR data
aircraft_data = xr.open_dataset(Obs_DIR + sfmr_file)

fp_time_min = 14000 #start time closest to first pass through storm center
fp_time_max = 14500 #end time closest to first pass through storm center
sp_time_min = 20000 #start time closest to second pass through storm center
sp_time_max = 21000 #end time closest to second pass through storm center

fp_wind_speed, fp_lat, fp_lon = find_min_wind_speed_location(aircraft_data, fp_time_min, fp_time_max)
sp_wind_speed, sp_lat, sp_lon = find_min_wind_speed_location(aircraft_data, sp_time_min, sp_time_max)

# Find the closest time for the given fp_lon and fp_lat
#AWO
awo_fp_closest_time = find_closest_time(awo_track, fp_lon, fp_lat)
awo_ws_fp_closest_time = find_closest_time(awo_ws_track, fp_lon, fp_lat)

# AWO_ws
awo_sp_closest_time = find_closest_time(awo_track, sp_lon, sp_lat)
awo_ws_sp_closest_time = find_closest_time(awo_ws_track, sp_lon, sp_lat)

print(awo_fp_closest_time, awo_sp_closest_time)
print(awo_ws_fp_closest_time, awo_ws_sp_closest_time)

# Get the file names for the closest times
#AWO
# awo_wave_files_fp = 'umwmout_awo_' + awo_fp_closest_time.strftime('%Y-%m-%d_%H:%M:%S') + '.nc'
# awo_wave_files_sp = 'umwmout_awo_' + awo_sp_closest_time.strftime('%Y-%m-%d_%H:%M:%S') + '.nc'
awo_wave_files_fp = 'wrfout_awo_d02_' + awo_fp_closest_time.strftime('%Y-%m-%d_%H:%M:%S') 
awo_wave_files_sp = 'wrfout_awo_d02_' + awo_sp_closest_time.strftime('%Y-%m-%d_%H:%M:%S')


# AWO_ws
# awo_ws_wave_files_fp = 'umwmout_awo_ws_' + awo_ws_fp_closest_time.strftime('%Y-%m-%d_%H:%M:%S') + '.nc'
# awo_ws_wave_files_sp = 'umwmout_awo_ws_' + awo_ws_sp_closest_time.strftime('%Y-%m-%d_%H:%M:%S') + '.nc'
awo_ws_wave_files_fp = 'wrfout_awo-ws_d02_' + awo_ws_fp_closest_time.strftime('%Y-%m-%d_%H:%M:%S')
awo_ws_wave_files_sp = 'wrfout_awo-ws_d02_' + awo_ws_sp_closest_time.strftime('%Y-%m-%d_%H:%M:%S')
#Accessing Data
#AWO
awo_wave_data_fp = xr.open_dataset(Model_DIR + awo_wave_files_fp)
awo_wave_data_sp = xr.open_dataset(Model_DIR + awo_wave_files_sp)
awo_wave_data_fp['wspd'] = np.sqrt(awo_wave_data_fp['U10']**2 + awo_wave_data_fp['V10']**2) #Computing wind speed from U10 and V10
awo_wave_data_sp['wspd'] = np.sqrt(awo_wave_data_sp['U10']**2 + awo_wave_data_sp['V10']**2) #Computing wind speed from U10 and V10

#AWO-ws
awo_ws_wave_files_fp = xr.open_dataset(Model_DIR + awo_ws_wave_files_fp)
awo_ws_wave_files_sp = xr.open_dataset(Model_DIR + awo_ws_wave_files_sp)
awo_ws_wave_files_fp['wspd'] = np.sqrt(awo_ws_wave_files_fp['U10']**2 + awo_ws_wave_files_fp['V10']**2) #Computing wind speed from U10 and V10
awo_ws_wave_files_sp['wspd'] = np.sqrt(awo_ws_wave_files_sp['U10']**2 + awo_ws_wave_files_sp['V10']**2) #Computing wind speed from U10 and V10

#Storm Info During First Pass
awo_storm_centers_fp = getStormCenter(plain2datetime(awo_fp_closest_time.strftime('%Y%m%d%H')),
                                      awo_track) #Identifying the track center at this particular time
#AWO
awo_storm_dir_fp = getStormDirection(plain2datetime(awo_fp_closest_time.strftime('%Y%m%d%H')), 
                                     awo_track) #Finding the track center
awo_x_wave_fp, awo_y_wave_fp = latlon2xyStormRelative(awo_wave_data_fp['XLONG'][0], 
                                                awo_wave_data_fp['XLAT'][0], 
                                                awo_storm_centers_fp[0], 
                                                awo_storm_centers_fp[1], 
                                                dir=awo_storm_dir_fp)


#AWO_ws
awo_ws_storm_centers_fp = getStormCenter(plain2datetime(awo_ws_fp_closest_time.strftime('%Y%m%d%H')),
                                         awo_ws_track) #Identifying the track center at this particular time
awo_ws_storm_dir_fp = getStormDirection(plain2datetime(awo_ws_fp_closest_time.strftime('%Y%m%d%H')), 
                                        awo_ws_track) #Finding the track center
awo_ws_x_wave_fp, awo_ws_y_wave_fp = latlon2xyStormRelative(awo_ws_wave_files_fp['XLONG'][0], 
                                                      awo_ws_wave_files_fp['XLAT'][0], 
                                                      awo_ws_storm_centers_fp[0], 
                                                      awo_ws_storm_centers_fp[1], 
                                                      dir=awo_ws_storm_dir_fp)

#Storm Info During Second Pass
awo_storm_centers_sp = getStormCenter(plain2datetime(awo_sp_closest_time.strftime('%Y%m%d%H')),
                                      awo_track) #Identifying the track center at this particular time
#AWO
awo_storm_dir_sp = getStormDirection(plain2datetime(awo_sp_closest_time.strftime('%Y%m%d%H')), 
                                     awo_track) #Finding the storm direction
awo_x_wave_sp, awo_y_wave_sp = latlon2xyStormRelative(awo_wave_data_sp['XLONG'][0], 
                                                awo_wave_data_sp['XLAT'][0], 
                                                awo_storm_centers_sp[0], 
                                                awo_storm_centers_sp[1], 
                                                dir=awo_storm_dir_sp)


#AWO_ws
awo_ws_storm_centers_sp = getStormCenter(plain2datetime(awo_ws_sp_closest_time.strftime('%Y%m%d%H')),
                                         awo_ws_track) #Identifying the track center at this particular time
awo_ws_storm_dir_sp = getStormDirection(plain2datetime(awo_ws_sp_closest_time.strftime('%Y%m%d%H')), 
                                        awo_ws_track) #Finding the storm direction
awo_ws_x_wave_sp, awo_ws_y_wave_sp = latlon2xyStormRelative(awo_ws_wave_files_sp['XLONG'][0], 
                                                      awo_ws_wave_files_sp['XLAT'][0], 
                                                      awo_ws_storm_centers_sp[0], 
                                                      awo_ws_storm_centers_sp[1], 
                                                      dir=awo_ws_storm_dir_sp)

#Obs
sfmr_x_fp, sfmr_y_fp = latlon2xyStormRelative(aircraft_data['LON'], aircraft_data['LAT'], 
                                        fp_lon, fp_lat, dir=awo_storm_dir_fp)

sfmr_x_sp, sfmr_y_sp = latlon2xyStormRelative(aircraft_data['LON'], aircraft_data['LAT'],
                                        sp_lon, sp_lat, dir=awo_storm_dir_sp)

ns_xmin = -5
ns_xmax = 5
ns_ymin = -200
ns_ymax = 200

ew_xmin = -200 
ew_xmax = 200
ew_ymin = -5
ew_ymax = 5

# AWO
awo_y_values_ns, awo_wspd_values_ns = extract_wspd_values_north_south(awo_x_wave_fp, awo_y_wave_fp, awo_wave_data_fp,
                                                ns_xmin, ns_xmax, ns_ymin, ns_ymax)

awo_x_values_ew, awo_wspd_values_ew = extract_wspd_values_east_west(awo_x_wave_sp, awo_y_wave_sp, awo_wave_data_sp,
                                                ew_xmin, ew_xmax, ew_ymin, ew_ymax)

#AWO_ws
awo_ws_y_values_ns, awo_ws_wspd_values_ns = extract_wspd_values_north_south(awo_ws_x_wave_fp, awo_ws_y_wave_fp, awo_ws_wave_files_fp,
                                                ns_xmin, ns_xmax, ns_ymin, ns_ymax)
awo_ws_x_values_ew, awo_ws_wspd_values_ew = extract_wspd_values_east_west(awo_ws_x_wave_fp, awo_ws_y_wave_fp, awo_ws_wave_files_fp,
                                                ew_xmin, ew_xmax, ew_ymin, ew_ymax)


x_dist = np.arange(-200, 201, 1)
y_dist = x_dist

interp_awo_wspd_smooth_ew = np.interp(x_dist, awo_x_values_ew, awo_wspd_values_ew)
interp_awo_ws_wspd_smooth_ew = np.interp(x_dist, awo_ws_x_values_ew, awo_ws_wspd_values_ew)

interp_awo_wspd_smooth_ns = np.interp(y_dist, awo_y_values_ns, awo_wspd_values_ns)
interp_awo_ws_wspd_smooth_ns = np.interp(y_dist, awo_ws_y_values_ns, awo_ws_wspd_values_ns)

# Apply Gaussian smoothing
awo_wspd_values_smooth_ns = gaussian_filter(interp_awo_wspd_smooth_ns, sigma=1)
awo_wspd_values_smooth_ew = gaussian_filter(interp_awo_wspd_smooth_ew, sigma=1)

awo_ws_wspd_values_smooth_ns = gaussian_filter(interp_awo_ws_wspd_smooth_ns, sigma=1)
awo_ws_wspd_values_smooth_ew = gaussian_filter(interp_awo_ws_wspd_smooth_ew, sigma=1)

#Masking East West winds so that we can get the North South winds
masked_sws_ns = np.where((sfmr_x_fp > -25) & (sfmr_x_fp < 25), aircraft_data['SWS'].values, np.nan)

masked_sws_ns = np.where((sfmr_x_fp > -50) & (sfmr_x_fp < 50) & (sfmr_y_fp > 20) & 
                         (sfmr_y_fp < 50) & (masked_sws_ns >= 0) & (masked_sws_ns <= 45), 
                         np.nan, masked_sws_ns)

#Masking North South winds so that we can get the East West winds
masked_sws_ew = np.where((sfmr_y_sp > -25) & (sfmr_y_sp < 25), aircraft_data['SWS'].values, 
                         np.nan)

masked_sws_ew = np.where((sfmr_x_sp > 0) & (sfmr_x_sp < 20) & (sfmr_y_sp > -25) & 
                        (sfmr_y_sp < 30) & (masked_sws_ew >= 25), 
                        np.nan, masked_sws_ew)

# wspd_diff_ew = awo_wspd_values_smooth_ew - awo_ws_wspd_values_smooth_ew
# wspd_diff_ns = awo_wspd_values_smooth_ns - awo_ws_wspd_values_smooth_ns


# interp_sws_ns = np.interp(y_dist, sfmr_y_fp, masked_sws_ns)
# print(interp_sws_ns)
# print(sfmr_x_fp[(sfmr_x_fp > -200) & (sfmr_x_fp < 200)])
# print(sfmr_y_fp)

#Identify values within -200 to 200 km in the East-West direction
mask_ew = (sfmr_x_fp >= -200) & (sfmr_x_fp <= 200)
sfmr_x_fp = sfmr_x_fp[mask_ew]
masked_sws_ew = masked_sws_ew[mask_ew]
# print(sfmr_x_fp.shape, masked_sws_ew.shape)

#Identify values within -200 to 200 km in the North-South direction
mask_ns = (sfmr_y_fp >= -200) & (sfmr_y_fp <= 200)
sfmr_y_fp = sfmr_y_fp[mask_ns]
masked_sws_ns = masked_sws_ns[mask_ns]

sws_ew_rol_avg = pd.Series(masked_sws_ew).rolling(window=49, min_periods=49, center=True).mean()
print(len(sws_ew_rol_avg))

plt.plot(sws_ew_rol_avg)
plt.show()

# interp_sws_ew = np.interp(x_dist, np.sort(sfmr_x_fp), masked_sws_ew)

print(len(x_dist))

# print(sfmr_x_fp.values.max(), sfmr_x_fp.values.min())

# wspd_obs_awo_ew = interp_sws_ew - awo_wspd_values_smooth_ew
# wspd_obs_awo_ws_ew = interp_sws_ew - awo_ws_wspd_values_smooth_ew
# print(wspd_obs_awo_ew)

xmin_sfmr = -200
xmax_sfmr = 200
ymin_sfmr = 0
ymax_sfmr = 20

xmin = -250
xmax = 250
ymin = xmin
ymax = xmax

ticks=np.arange(-250, 300, 50)

xloc, yloc = 0.05, 0.88
#Figure Settings
fontsize=6
labelsize=6
gridsize =(1,2)
s=10
lw=1
# zorder=200
# size=50

#Colorbar Settings
wspd_levels = np.arange(5,57.5,2.5)
CMAPS = cmaps.MPL_gist_rainbow_r
shrink=0.70
aspect=18
width=0.5
length=4

labelpad = 0.25

#Plotting the data
fig = plt.figure(figsize=(4, 4))

#AWO
ax3 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)

# ax3.plot(x_dist, wspd_obs_awo_ew, c='k', lw=lw, label='CTL')
# ax3.plot(x_dist, awo_ws_wspd_values_smooth_ew, c='blue', lw=lw, label='EXP')
ax3.plot(sfmr_x_fp, interp_sws_ew, c='k', lw=lw, label='SFMR')



# ax3.set_xlim(xmin_sfmr, xmax_sfmr)
# ax3.set_ylim(ymin_sfmr, ymax_sfmr) 

ax3.set_xlabel('Distance to Storm Center (km)', fontsize=fontsize)
ax3.set_ylabel('Wind Speed (m/s)', fontsize=fontsize)
ax3.set_title('West - East', fontsize=fontsize)
ax3.grid(True, linestyle='--', linewidth=0.5)
# ax3.set_yticks(np.arange(-20, 25, 10))
# ax3.set_yticklabels(np.arange(-20, 25, 10), fontsize=fontsize)
# ax3.set_xticks(np.arange(-200, 250, 100))
# ax3.set_xticklabels(np.arange(-200, 250, 100), fontsize=fontsize)
add_corner_label(ax3, xloc, yloc, '(c)', fontsize)
ax3.legend(loc='upper right', fontsize=4)

ax4 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)

# ax4.plot(y_dist, wspd_diff_ns, lw=lw, c='k', label='CTL')
# ax4.plot(y_dist, awo_ws_wspd_values_smooth_ns, lw=lw, c='blue', label='EXP')
ax4.plot(sfmr_y_fp, masked_sws_ns, lw=lw, c='k', label='SFMR')

# ax4.set_xlim(xmin_sfmr, xmax_sfmr)
# ax4.set_ylim(ymin_sfmr, ymax_sfmr)

ax4.set_xlabel('Distance to Storm Center (km)', fontsize=fontsize)
ax4.set_ylabel('Wind Speed (m/s)', fontsize=fontsize)
ax4.set_title('North - South', fontsize=fontsize)
# ax4.set_yticks(np.arange(-20, 25, 10))
# ax4.set_yticklabels(np.arange(-20, 25, 10), fontsize=fontsize)
# ax4.set_xticks(np.arange(-200, 250, 100))
# ax4.set_xticklabels(np.arange(-200, 250, 100), fontsize=fontsize)
ax4.grid(True, linestyle='--', linewidth=0.5)
add_corner_label(ax4, xloc, yloc, '(d)', fontsize)
fig.tight_layout(pad=0, w_pad=0, h_pad=0)


plt.savefig(PNG + 'asymmetry_diff.png', dpi=300, bbox_inches='tight', 
            facecolor='w', edgecolor='w', transparent=False)

# plt.show()