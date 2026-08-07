import numpy as np
import pandas as pd
import xarray as xr
import glob
import os
import datetime 
from datetime import datetime,timedelta

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cmaps

from helpers import *

SIM_PATH = '/home/disk/orca/adaley17/Research/Stress_Separation/Hurricane_Earl/Data/awo5-ws_2010082700_gfs_3.7.1/'
PNG = '/home/disk/orca/adaley17/orca3/adaley17/Projects/Air_Sea_Momentun_Exchange_Open_Ocean/Figures/'

earl_wrf_dom1 = 'wrfout_d01_2010-09-01_12:00:00' # WRF output file for Domain 1 at the specified time
earl_wrf_dom2 = 'wrfout_d02_2010-09-01_12:00:00' # WRF output file for Domain 2 at the specified time
earl_wrf_dom3 = 'wrfout_d03_2010-09-01_12:00:00' # WRF output file for Domain 3 at the specified time
earl_hyc_dom = 'archv.2010_243_19.nc' # HYCOM output file for the specified time (2010-09-01 12:00:00)
earl_dwm_dom = 'umwmout_2010-09-01_12:00:00.nc' # DWM output file for the specified time (2010-09-01 12:00:00)

earl_wrf_01 = xr.open_dataset(SIM_PATH + earl_wrf_dom1) # WRF Domain 1
earl_wrf_02 = xr.open_dataset(SIM_PATH + earl_wrf_dom2) # WRF Domain 2
earl_wrf_03 = xr.open_dataset(SIM_PATH + earl_wrf_dom3) # WRF Domain 3
earl_hyc = xr.open_dataset(SIM_PATH + earl_hyc_dom) # HYCOM Domain
earl_dwm = xr.open_dataset(SIM_PATH + 'output/' + earl_dwm_dom) # DWM Domain

earl_wrf_01['WSPD'] = np.sqrt(earl_wrf_01['U10']**2 + earl_wrf_01['V10']**2) #Computing Wind Speed from U and V components


earl_hycom_min_lon = float(np.nanmin(earl_hyc['longitude'])) #Beryl hycom domain min lon
earl_hycom_max_lon = float(np.nanmax(earl_hyc['longitude'])) #Beryl hycom domain max lon

earl_hycom_min_lat = float(np.nanmin(earl_hyc['latitude'])) #Beryl hycom domain min lat
earl_hycom_max_lat = float(np.nanmax(earl_hyc['latitude'])) #Beryl hycom domain max lat

lw=1

#Colorbar Settings
shrink=0.40#Colorbar Shrink
aspect=15 #Colorbar Aspect
width=0.5 #Colorbar Width
length=4 #Colorbar Length
labelsize=6 #Colorbar Label Size

#Cartopy Settings
min_lat = 4
max_lat = 45.5
min_lon = -106
max_lon = -20
CMAPS=cmaps.MPL_jet_r
proj = ccrs.PlateCarree()
# proj._threshold /= 20.  #allows fine grain plot
coast = cfeature.GSHHSFeature(scale='high', levels=[1,], edgecolor='k')
lakes = cfeature.GSHHSFeature(scale='high', levels=[2,], edgecolor='face',facecolor='grey')
plot_area = [min_lon,max_lon, min_lat, max_lat]


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1, projection=proj)

wspd = ax.contourf(earl_wrf_01['XLONG'][0], earl_wrf_01['XLAT'][0], earl_wrf_01['WSPD'][0],
            levels=np.arange(0, 37.5, 0.5), cmap=CMAPS, extend='max',
            transform=ccrs.PlateCarree())

hcb = fig.colorbar(wspd, shrink=shrink, aspect=aspect, ax=ax, pad=0.02)
hcb.set_label('Wind Speed (m/s)', fontsize=labelsize)
hcb.ax.tick_params(color='k', length=length, width=width, labelsize=labelsize, pad=0.002)

Cartopy_Features(ax, 12, plot_area, 4, 4, 'k') #Add Cartopy Features

# Domain01
plt.plot(earl_wrf_01['XLONG'][0][:,0], earl_wrf_01['XLAT'][0][:,0], '-', color='k',  lw=lw, transform=ccrs.PlateCarree(), label='WRF Dom1')
plt.plot(earl_wrf_01['XLONG'][0][:,-1],  earl_wrf_01['XLAT'][0][:,-1], '-', color='k',  lw=lw, transform=ccrs.PlateCarree(), )
plt.plot(earl_wrf_01['XLONG'][0][0,:],  earl_wrf_01['XLAT'][0][0,:], '-', color='k',  lw=lw, transform=ccrs.PlateCarree(), )
plt.plot(earl_wrf_01['XLONG'][0][-1,:],  earl_wrf_01['XLAT'][0][-1,:], '-', color='k',  lw=lw, transform=ccrs.PlateCarree(), )


#Domain02
plt.plot(earl_wrf_02['XLONG'][0][:,0],  earl_wrf_02['XLAT'][0][:,0], '--', color='k',  lw=lw, transform=ccrs.PlateCarree(), label='WRF Dom2')
plt.plot(earl_wrf_02['XLONG'][0][:,-1],  earl_wrf_02['XLAT'][0][:,-1], '--', color='k',  lw=lw, transform=ccrs.PlateCarree())
plt.plot(earl_wrf_02['XLONG'][0][0,:],  earl_wrf_02['XLAT'][0][0,:], '--', color='k',  lw=lw, transform=ccrs.PlateCarree())
plt.plot(earl_wrf_02['XLONG'][0][-1,:],  earl_wrf_02['XLAT'][0][-1,:], '--', color='k',  lw=lw, transform=ccrs.PlateCarree())


# Domain03
plt.plot(earl_wrf_03['XLONG'][0][:,0],  earl_wrf_03['XLAT'][0][:,0], '-.', color='k',  transform=ccrs.PlateCarree(), label='WRF DOM3')
plt.plot(earl_wrf_03['XLONG'][0][:,-1],  earl_wrf_03['XLAT'][0][:,-1], '-.', color='k',  transform=ccrs.PlateCarree() )
plt.plot(earl_wrf_03['XLONG'][0][0,:],  earl_wrf_03['XLAT'][0][0,:], '-.', color='k',  transform=ccrs.PlateCarree() )
plt.plot(earl_wrf_03['XLONG'][0][-1,:],  earl_wrf_03['XLAT'][0][-1,:], '-.', color='k',  transform=ccrs.PlateCarree())

# #DWM Domain
# #DWM Domain
plt.plot(earl_dwm['lon'][0][:,0],  earl_dwm['lat'][0][:,0], '--', color='cyan',  lw=lw, transform=ccrs.PlateCarree(), label='DWM DOM')
plt.plot(earl_dwm['lon'][0][:,-1],  earl_dwm['lat'][0][:,-1], '--', color='cyan',  lw=lw, transform=ccrs.PlateCarree())
plt.plot(earl_dwm['lon'][0][0,:],  earl_dwm['lat'][0][0,:], '--', color='cyan',  lw=lw, transform=ccrs.PlateCarree())
plt.plot(earl_dwm['lon'][0][-1,:],  earl_dwm['lat'][0][-1,:], '--', color='cyan',  lw=lw, transform=ccrs.PlateCarree())


# #Hycom Domain
plt.plot([earl_hycom_min_lon, earl_hycom_min_lon], [earl_hycom_min_lat, earl_hycom_max_lat], ':', color='k',  lw=lw, transform=ccrs.PlateCarree(), label='HYC DOM')
plt.plot([earl_hycom_min_lon, earl_hycom_max_lon], [earl_hycom_max_lat, earl_hycom_max_lat], ':', color='k',  lw=lw, transform=ccrs.PlateCarree())
plt.plot([earl_hycom_max_lon, earl_hycom_max_lon], [earl_hycom_min_lat, earl_hycom_max_lat], ':', color='k',  lw=lw, transform=ccrs.PlateCarree())
plt.plot([earl_hycom_min_lon, earl_hycom_max_lon], [earl_hycom_min_lat, earl_hycom_min_lat], ':', color='k',  lw=lw, transform=ccrs.PlateCarree())

ax.legend(loc='upper right', fontsize=4, shadow=True)

plt.savefig(PNG + 'fig1_earl_domain_ESS.png', dpi=400, bbox_inches='tight',
                 facecolor='w', transparent=False)

print("We are Finish")