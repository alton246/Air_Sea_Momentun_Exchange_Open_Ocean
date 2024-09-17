import numpy as np
import xarray as xr
import datetime as dt
import cftime
import matplotlib.pyplot as plt
import cmaps
import math
import numpy.ma as ma

from helpers import *

#We first need to get the paths where the data is located
PATH = '/home/disk/orca/adaley17/Research/Stress_Separation/Hurricane_Earl/Data/06HR_Composite_Averages_QC/new_all/'
PNG = '/home/disk/orca/adaley17/my_stuff/Publications/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'

#We also need the file name
file = 'storm_relative_qualit_controlled_data.*.nc'

#We need to access the data 
comp_data = xr.open_mfdataset(PATH + file)

#We need to extract from the data the period we are interested in 
data_subset = comp_data.sel(time=slice('2010-08-30T12:00:00.000000000', '2010-09-01T06:00:00.000000000'))

#Dividing MLD by 9806 to change to meters
data_subset['mld_awo_qc'] = data_subset['mld_awo_qc']/9806
data_subset['mld_awo_ws_qc'] = data_subset['mld_awo_ws_qc']/9806


# Creating a mask to identify region of cold wake

sst_array = data_subset['sst_awo_ws_qc'][7].to_masked_array()  #Changing data to a masked array

sst_criteria = (data_subset['sst_awo_ws_qc'][7].values > 28.4) #Criteria to identify the data outside the points of interest

sst_mask = ma.masked_where(sst_criteria, sst_array) #Executing criteria on data to create the masked data

sst_masked_ind = np.where(sst_mask.mask) #Indentifying the index of the data to be excluded

sst_mask[sst_masked_ind] = np.nan #Converting the points outside the cold wake to nans

awo_ws_sst_masked_qc = sst_mask.data # Changing masked data to regular data format

data1 = data_subset['swh_awo_qc']
data2 = data_subset['swh_awo_ws_qc']
n = len(data_subset['time'])
x = np.arange(1,n+1,1)

#Etracting Waves ahead and behind of the storm
#SWH
swh_ahead, swh_medians_ahead = extract_data_ahead_of_storm(data1, data2, n)
swh_behind, swh_medians_behind = extract_data_behind_of_storm(data1, data2, n)

#Computing the best fit line
swh_ahead_bf_line, swh_ahead_slope, swh_ahead_intercept = generate_best_fit_line(x, swh_medians_ahead, 1)
swh_behind_bf_line, swh_behind_slope, swh_behind_intercept = generate_best_fit_line(x, swh_medians_behind, 1)

#Extracting Data SST and SFC data within the cold wake.
#SST
sst_data, sst_medians = extract_masked_data_of_interest(data_subset['sst_awo_qc'], data_subset['sst_awo_ws_qc'], 
                                                      sst_masked_ind, n)

#Computing the best fit line
sst_bf_line, sst_slope, sst_intercept = generate_best_fit_line(x, sst_medians, 1)


#SFC
sfc_data, sfc_medians = extract_masked_data_of_interest(data_subset['sfc_awo_qc'], data_subset['sfc_awo_ws_qc'], 
                                                      sst_masked_ind, n)

#Computing the best fit line
sfc_bf_line, sfc_slope, sfc_intercept = generate_best_fit_line(x, sfc_medians, 1)

awo_data = ['swh_awo_qc', 'sfc_awo_qc', 'sst_awo_qc']
awo_ws_data = ['swh_awo_ws_qc', 'sfc_awo_ws_qc', 'sst_awo_ws_qc']
time_data = data_subset['time'].values
n = len(data_subset['time'])
x = np.arange(1,n+1,1)
num_var = len(awo_data)

#Color Options
color_options = ['red', 'gray', 'green']

#Figure Options
ylabels = ['$H_{s}$ Diff (m)', '$U$ Diff $(m/s)$', '$SST$ Diff $(^{\circ}C)$ ']
ylim_max = [1.0, 1.0, 0.5, 0.2, 10]
ylim_min = [-1.0, -1.0, -0.5, -0.2, -10]
median_pos = [0.80, 0.15, 0.45, 0.45, 9.0]


#Time of interest
time_info = pd.date_range('2010-08-30 15:00:00',
                          '2010-09-01 09:00:00',
                          freq='3H')
# Converting time
formatted_dates = [pd.Timestamp(date).strftime('%m-%d %H') for date in time_info[::2]]




#Creating Plot
edge_color='k'
lw = 1
alpha=0.8
fontsize=8
labelsize=8
rotation=45
zorder=200
txt_color='magenta'
length=1.5
width=1
labelpad=0.005
markersize=10
swh_ymin=-0.5
swh_ymax=0.5
sfc_ymin=-0.25
sfc_ymax=0.25
fontsize_title=10
title_pad=0.05
y_pos = 0.85
x_pos = 0.93
cl_fontsize=6

gridsize = (3, 1)
fig = plt.figure(figsize=(4, 5))


# AWO 0000 UTC
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
# ax1_twin = ax1.twinx()

#Plots
ax1.scatter(x, swh_medians_ahead, color=color_options[0], s=markersize, zorder=100)
ax1.axhline(y=0, color='blue', linestyle='--', linewidth=lw, label='y=0') #Horizontal line
ax1.plot(x, swh_ahead_bf_line, color=color_options[0], lw=lw, zorder=100)


#Bulid a function
ax1.text(x[0], -0.4, str(np.round(np.nanmean(swh_medians_ahead),5)) + ' m', 
             fontsize=fontsize, zorder=zorder, color=color_options[0])
ax1.set_xticks(x)
ax1.set_xticklabels(formatted_dates)

#Axis labels
add_axis_label_to_single_axis_summary_plot(ax1, ylabels[0], fontsize, labelpad, 
                                           labelsize, rotation, swh_ymin, swh_ymax, length, width)

ax1.set_title('Ahead', loc='left', fontsize=fontsize_title, fontweight='bold', pad=title_pad)

add_corner_label(ax1, x_pos, y_pos, '(a)', fontsize=cl_fontsize)


ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)

#Plots
ax2.scatter(x, swh_medians_behind, color=color_options[0], s=markersize, zorder=100)
# ax2_twin.scatter(x, wspd_medians_FR, color=color_options[1], s=markersize, zorder=100)
ax2.axhline(y=0, color='blue', linestyle='--', linewidth=lw, label='y=0') #Horizontal line
ax2.plot(x, swh_behind_bf_line, color=color_options[0], lw=lw, zorder=100)
# ax2_twin.plot(x, wspd_FR_bf_line, color=color_options[1], lw=lw, zorder=100, 
#          label=f'Best-fit line: y = {wspd_FR_slope:.2f}x + {wspd_FR_intercept:.2f}')

#Bulid a function
ax2.text(x[0], -0.4, str(np.round(np.nanmean(swh_medians_behind),5)) + ' m', 
             fontsize=fontsize, zorder=zorder, color=color_options[0])


#Format time for x labels
ax2.set_xticks(x)
ax2.set_xticklabels(formatted_dates)

#Axis labels
add_axis_label_to_single_axis_summary_plot(ax2, ylabels[0], fontsize, labelpad, 
                                           labelsize, rotation, swh_ymin, swh_ymax, length, width)
ax2.set_title('Behind', loc='left', fontsize=fontsize_title, fontweight='bold', pad=title_pad)

add_corner_label(ax2, x_pos, y_pos, '(b)', fontsize=cl_fontsize)

#Build a function
#Legend
# add_legend_to_summary_plot(ax2, color_options[0], color_options[1], '$H_{s}$', '$U_{10}$', 'lower left', 4.5)


ax3 = plt.subplot2grid(gridsize, (2, 0), colspan=1, rowspan=1)
ax3_twin = ax3.twinx()

#Plots
ax3.scatter(x, sfc_medians, color=color_options[1], s=markersize, zorder=100)
ax3_twin.scatter(x, sst_medians, color=color_options[2], s=markersize, zorder=100)
plt.axhline(y=0, color='blue', linestyle='--', linewidth=lw, label='y=0') #Horizontal line
ax3.plot(x, sfc_bf_line, color=color_options[1], lw=lw, zorder=100, 
         label=f'Best-fit line: y = {sfc_slope:.2f}x + {sfc_intercept:.2f}')
ax3_twin.plot(x, sst_bf_line, color=color_options[2], lw=lw, zorder=100, 
         )

ax3.text(x[0], -0.18, str(np.round(np.nanmean(sfc_medians),5)) + ' m/s', 
             fontsize=fontsize, zorder=zorder, color=color_options[1])
ax3.text(x[0], -0.24, str(np.round(np.nanmean(sst_medians),5)) + '$^{\circ}C$', 
             fontsize=fontsize, zorder=zorder, color=color_options[2])

#Format time for x labels
ax3.set_xticks(x)
ax3.set_xticklabels(formatted_dates)

#Axis labels
add_axis_label_to_dub_axis_summary_plot(ax3, ax3_twin, ylabels[1], ylabels[2], fontsize, labelpad, 
                               labelsize, color_options[1], color_options[2], rotation, sfc_ymin, sfc_ymax, length, width)
ax3.set_title('Cold Wake', loc='left', fontsize=fontsize_title, fontweight='bold', pad=title_pad)

add_corner_label(ax3, x_pos, y_pos, '(c)', fontsize=cl_fontsize)

#Legend
add_legend_to_summary_plot(ax3, color_options[1], color_options[2], 'U', 'SST', 'upper left', 4.5)


fig.tight_layout(pad=0, w_pad=0, h_pad=0)

plt.savefig(PNG + 'fig5_summary_swh_sst_sfc_contour.png', dpi=300, bbox_inches='tight',
                facecolor='w', transparent=False)


print('The script was executed succesfully!!')
