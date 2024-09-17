import numpy as np
import pandas as pd
import datetime
import csv
import matplotlib.pyplot as plt
import math

from matplotlib.dates import num2date,date2num
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
from datetime import datetime,timedelta
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import xarray as xr


class atcf_csv():
    """
    ATCF class that takes data from .csv file.
    """
    def __init__(self,path):
        self.time = []
        self.lat  = []
        self.lon  = []
        self.wspd = []
        self.mslp = []
        # Read data
        df = pd.read_csv(path,parse_dates=['date_and_time'])
        self.time = np.array([t.to_pydatetime() for t in df['date_and_time']])
        self.lat  = df['lat'].to_numpy()    # [deg N]
        self.lon  = df['lon'].to_numpy()    # [deg E]
        self.wspd = df['max_wind_kts'].to_numpy()*0.514444    # [m s-1]
        self.mslp = df['min_slp_mb'].to_numpy()    # [mb]

class extract_ibt_data():
    def __init__(self, num_entries, data):

        j = 0 # Looping Variable used to fill the

        #Empty arrays to store respective data from hurdat file
        self.lats = np.zeros((num_entries,1), dtype=float)
        self.lons = np.zeros((num_entries,1), dtype=float)
        self.colors = np.zeros((num_entries,1), dtype=object)
        self.max_ws = np.zeros((num_entries,1), dtype=float)
        self.max_pres = np.zeros((num_entries,1), dtype=float)
        self.date = np.zeros((num_entries,1), dtype=datetime)


        for rows in range(len(data)):
            if data[rows][5] == 'EARL' and data[rows][1] == '2010':
            
                # Extracting dates from the NHC best track
                self.date[j,:] = datetime.strptime(data[rows][6], '%Y-%m-%d %H:%M:%S')

                # Extracting latitudes and longitudes of storm center 
                self.lats[j,:] = data[rows][8]
                self.lons[j,:] = data[rows][9]

                # Extracting maximum wind speeds 
                self.max_ws[j,:] = float(data[rows][23]) / 1.944

                # Assinging colors to markers based on the intensity of the storm
                self.colors[j,:] = ReturnColorIntensity(float(data[rows][23]))

                # Extracting minimum sea level pressure 
                self.max_pres[j,:] = float(data[rows][24])

                j = j + 1

def array_scrub(array):
    scrubbed_array = []
    for e in array:
        value = e.strip()
        if value == '-999':
            scrubbed_array += [""]
        elif value == 'NaN':
            scrubbed_array += [""]
        else:
            scrubbed_array += [value]
    return scrubbed_array

def make_wind_radii_obj(input_array):
    radii_array = array_scrub(input_array)
    # initialize empty object for wind radii extents
    wind_radii_obj = {}
    header = ['34kt', '50kt', '64kt']

    # quadrants follow NE, SE, SW, NW convention
    i = 0
    for head in header:
        wind_radii_obj[head] = make_wind_coords_obj(radii_array[i : i + 4])
        i += 4
    return wind_radii_obj


def make_wind_coords_obj(input_array):
    header = ['NE', 'SE', 'SW', 'NW']
    wind_coords = {}
    i = 0
    for e in header:
        wind_coords[e] = input_array[i]
        i += 1
    return wind_coords

def make_coords_obj(coords):
    coords_obj = {}
    coords_obj['lat'], coords_obj['lon'] = coords[0], coords[1]
    return coords_obj

def make_date_obj(date):
    date_obj = {}
    date_obj['year'] = int(date[0:4])
    date_obj['month'] = int(date[4:6])
    date_obj['day'] = int(date[-2:])
    return date_obj

def plain2datetime(timeStr):
    """
    Returns a datetime instance based on 
    a time string in form YYYYMMDDhh.
    """
    YYYY = timeStr[0:4]
    # print(YYYY)
    MM = timeStr[4:6]
    DD = timeStr[6:8]
    hh = timeStr[8:10]
    return datetime(int(YYYY),int(MM),int(DD),int(hh),0,0)

def convert_latlon_to_decimal(coord):
    # Change coordinate notation to decimal
    # -N, E hemispheres are positive
    # -S, W hemispheres are negative
    
    if coord[-1] == 'N' or coord[-1] == 'E':
        lat_lon = float(coord[:-1])
    else:
        lat_lon = -1 * float(coord[:-1])               
    return lat_lon


def argminDatetime(time0,time):
    """
    Returns the index of datetime array time for which
    the datetime time0 is nearest.
    """
    return np.argmin(np.abs(date2num(time0)-date2num(time)))

def getStormDirection(time,track):
    """
    Given datetime and atcf track objects, this function returns
    storm propagation direction.
    """
    nup = argminDatetime(time+timedelta(hours=3),track.time)
    ndn = argminDatetime(time-timedelta(hours=3),track.time)

    dlon = track.lon[nup]-track.lon[ndn]
    dlat = track.lat[nup]-track.lat[ndn]

    return np.arctan2(dlat,dlon)

def getStormCenter(time,track):
    """
    Returns the (lon,lat) storm center fiven an atcf track instance.
    """
    ind = np.argmin(np.abs(date2num(time)-date2num(track.time)))
    return track.lon[ind],track.lat[ind]


def latlon2xyStormRelative(lon,lat,lon0,lat0,dir=0.5*np.pi):
    """
    Given (lon,lat) fields, returns (x,y) distance in meters 
    from (lon0,lat0). If dir [radians] is specified, rotates (x,y)
    so that positive y-axis is pointing in direction dir.
    """
    Re   = 6.371E6
    d2km = np.pi*Re/180.*1E-3
    d2r  = np.pi/180.
    Y = (lat-lat0)*d2km
    X = (lon-lon0)*d2km*np.cos(lat*d2r)
    x = X*np.sin(dir)-Y*np.cos(dir)
    y = X*np.cos(dir)+Y*np.sin(dir)
    return x,y

def extract_Storm_centers_over_period_time(date_of_int, track_data):
    storm_centers = np.zeros((len(date_of_int),2), dtype=float)
    for i in range(len(date_of_int)):
        lats_lons_awo = getStormCenter(date_of_int[i], track_data)
        storm_centers[i,0] = lats_lons_awo[0]
        storm_centers[i,1] = lats_lons_awo[1]

    return storm_centers

def Rotate_U_V_Vector_Comp(u, v, storm_dir):
    u_rot = u*np.cos(storm_dir-(np.pi/2)) + v*np.sin(storm_dir-(np.pi/2))
    v_rot = -u*np.sin(storm_dir-(np.pi/2)) + v*np.cos(storm_dir-(np.pi/2)) 

    return u_rot, v_rot
def counter(PATH, track_file, Storm_Name, Year):
    #File Reader
    with open(PATH+ track_file) as file:
        csvreader = csv.reader(file)
        new_rows = list(csvreader)

    # Entry Counter
    counter = 0
    for i in range(len(new_rows)):
        if new_rows[i][5] == Storm_Name and new_rows[i][1] == Year:
            counter = counter + 1

    return counter, new_rows

def ReturnColorIntensity(wspds):
    col = []
    wspd_kts = wspds
    # for i in range(len(wspd_kts)):
        # print(i, wspd_kts[i])
    if wspd_kts <= 33.0:
            # print(bt_data[rows][23])
        col.append('green')
        # Tropical storm
    elif wspd_kts >= 33.0 and wspd_kts <=63.0:
        # print(bt_data[rows][23])
        col.append('yellow')
    #Hurricane Category 1
    elif wspd_kts >= 64.0 and wspd_kts <83.0:
        # print(bt_data[rows][23])
        col.append('gold')
    # Hurricane category 2
    elif wspd_kts >= 83.0 and wspd_kts <96.0:
        col.append('orange')
    # Hurricane category 3
    elif wspd_kts >= 96.0 and wspd_kts <113.0:
        col.append('tomato')
    # Hurricane category 4
    elif wspd_kts >= 113.0 and wspd_kts <137.0:
        col.append('red')
    # Hurricane category 5
    elif wspd_kts >= 137.0:
        col.append('darkred')
    else:
        col.append(np.nan)

    if wspd_kts == ' ':
            col.append(np.nan)
    return col

def Cartopy_Features(axis, plot_area, xlocator, ylocator):
        """Plot area is a list in the form 
        [max_lon_track + 360.0, min_lon_track + 360.0, min_lat_track,max_lat_track]
        """
        #Setting coast, lakes and projection
        coast = cfeature.GSHHSFeature(scale='high', levels=[1,], edgecolor='k')
        lakes = cfeature.GSHHSFeature(scale='high', levels=[2,], edgecolor='face',facecolor='grey')
        axis.add_feature(cfeature.LAND, color='white')
        axis.add_feature(coast)
        axis.add_feature(lakes)
        axis.add_feature(cfeature.BORDERS, linewidth=0.75)
        # ax1.add_feature(cfeature.OCEAN)
        states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='110m')
        axis.add_feature(states_provinces, linewidth=0.75)
        axis.set_extent(plot_area, crs=ccrs.PlateCarree())

        gl = axis.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=0.5, color='k', alpha=0.7, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        # gl.left_labels = False
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 180,xlocator))
        gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, ylocator))
        gl.xlabel_style = {'size': 15, 'color': 'black'}
        gl.ylabel_style = {'size': 15, 'color': 'black'}

        return

def Cartopy_Features(axis, fontsize, plot_area, xlocator, ylocator):
        """Plot area is a list in the form 
        [max_lon_track + 360.0, min_lon_track + 360.0, min_lat_track,max_lat_track]
        """

def julian_to_regular(julian_day):
    # Calculate the base date for Julian day
    base_date = datetime(2010, 1, 1)
    # Add the number of days to the base date
    target_date = base_date + timedelta(days=julian_day)
    return target_date.strftime("%Y-%m-%d")

def add_storm_relative_essentials(axis):
    # axis.axvline(0, color='k', linestyle='solid', linewidth=0.5)
    # axis.axhline(0, color='k', linestyle='solid', linewidth=0.5)
    # axis.arrow(0, 0, 0, 250, width=1.5, color='k')

    axis.axvline(0, color='k', linestyle='solid', linewidth=0.5)
    axis.axhline(0, color='k', linestyle='solid', linewidth=0.5)
    axis.arrow(0, 0, 0, 175, width=0.75, color='k')
    first_circle_wspd = plt.Circle( (0, 0), 50, lw=0.5, fill = False, color='k')
    second_circle_wspd = plt.Circle( (0, 0), 150, lw=0.5, fill = False, color='k')
    third_circle_wspd = plt.Circle( (0, 0), 250, lw=0.5, fill = False, color='k')
    axis.add_artist(first_circle_wspd)
    axis.add_artist(second_circle_wspd)
    axis.add_artist(third_circle_wspd)

def getSWHSymmetry(dist_x_vals, dist_y_vals, data, min_x_dist, max_x_dist, min_y_dist, max_y_dist):
    dist_x_vals = dist_x_vals.values.flatten() #flattening data
    dist_y_vals = dist_y_vals.values.flatten() #flattening data

    #obtaining indices of distances between 0 and 100 km
    indices =np.where((dist_y_vals>=min_y_dist) & (dist_y_vals<=max_y_dist ) & 
         (dist_x_vals >= min_x_dist) & (dist_x_vals <= max_x_dist))
    
    dist_subset=dist_x_vals[indices]

    #Locating swh values that are between 0 and 100 km
    swh_vals = data.values.flatten() #flattening wave height data
    swh_subset= swh_vals[indices]

    return swh_subset, dist_subset

def EquatLenOfList(val1, val2):
    num_nans = np.abs(len(val1) - len(val2))

    if len(val1) < len(val2):
        new_val1 = val1.tolist() + [float('nan')] * num_nans
        new_val2 = val2

    elif len(val2) < len(val1):
        new_val1 = val2.tolist() + [float('nan')] * num_nans
        new_val2= val1

    else:
        new_val1 = val1.copy()
        new_val2 = val2.copy()

    return new_val1, new_val2 ;

def AlongTrackCrossSection(awo_x_dist, awo_y_dist, awo_temp, awo_ws_x_dist, awo_ws_y_dist, awo_ws_temp):

    awo_surf_temp_symmetry, awo_x_dist_subset_temp = getSWHSymmetry(awo_x_dist, awo_y_dist, awo_temp,
                                                            -250, 250, -1, 2)
    awo_ws_surf_temp_symmetry, awo_ws_x_dist_subset_temp = getSWHSymmetry(awo_ws_x_dist, awo_ws_y_dist, 
                                                                awo_ws_temp, -250, 250, -1, 2)

    awo_temp_tested, awo_ws_temp_tested = EquatLenOfList(awo_surf_temp_symmetry, awo_ws_surf_temp_symmetry)
    awo_dist_tested, awo_ws_dist_tested = EquatLenOfList(awo_x_dist_subset_temp, awo_ws_x_dist_subset_temp)

    return awo_temp_tested, awo_dist_tested, awo_ws_temp_tested, awo_ws_dist_tested

def GetBathyData(PATH, file):
    awo_ocn_data = xr.open_dataset(PATH + 'awo5_2010082700_gfs_3.7.1/' + file)

    return awo_ocn_data['bathymetry']

def Cartopy_Features(axis, fontsize, plot_area, xlocator, ylocator, linecolor):
        """Plot area is a list in the form 
        [max_lon_track + 360.0, min_lon_track + 360.0, min_lat_track,max_lat_track]
        """

        # pylab.rcParams['xtick.major.pad']='0'
        # pylab.rcParams['ytick.major.pad']='0'


        #Setting coast, lakes and projection
        coast = cfeature.GSHHSFeature(scale='high', levels=[1,], edgecolor='k')
        lakes = cfeature.GSHHSFeature(scale='high', levels=[2,], edgecolor='face',facecolor='grey')
        axis.add_feature(cfeature.LAND, color='white')
        axis.add_feature(coast)
        axis.add_feature(lakes)
        axis.add_feature(cfeature.BORDERS, linewidth=0.75)
        axis.add_feature(cfeature.BORDERS, linewidth=0.25)
        # ax1.add_feature(cfeature.OCEAN)
        states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='110m')
        axis.add_feature(states_provinces, linewidth=0.25)
        axis.set_extent(plot_area, crs=ccrs.PlateCarree())

        gl = axis.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=0.2, color=linecolor, alpha=0.7, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        # gl.bottom_labels = False
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 180,xlocator))
        gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, ylocator))
        gl.xlabel_style = {'size': fontsize, 'color': 'black'}
        gl.ylabel_style = {'size': fontsize, 'color': 'black'}
        gl.xpadding = 0.5
        gl.ypadding = 0.5

def Track_Legend(axis, fontsize, markersize, xlocator, ylocator, loc):
        #Handles for legend
        ts = mlines.Line2D([], [], color='yellow', marker='o', markersize=markersize, ls='', label='tropical strom')
        hur_cat1 = mlines.Line2D([], [], color='gold', marker='o', markersize=markersize, ls='', label='hurricane cat1')
        hur_cat2 = mlines.Line2D([], [], color='orange', marker='o', markersize=markersize, ls='', label='hurricane cat2')
        hur_cat3 = mlines.Line2D([], [], color='tomato', marker='o', markersize=markersize, ls='', label='hurricane cat3')
        hur_cat4 = mlines.Line2D([], [], color='red', marker='o', markersize=markersize, ls='',  label='hurricane cat4')
        # hur_cat5 = mlines.Line2D([], [], color='darkred', marker='o', ls='', label='hurricane cat5')
        bt = mlines.Line2D([], [], color='black', marker=' ', markersize=markersize, ls='solid', label='best track')
        awo_ws = mlines.Line2D([], [], color='cyan', marker='*', markersize=markersize,ls='solid', label='$AWO_{ws}$-$EXP$')
        awo = mlines.Line2D([], [], color='red', marker='s', markersize=markersize, ls='solid', label='$AWO$-$CTL$')
        ao = mlines.Line2D([], [], color='magenta', marker='s', markersize=markersize, ls='solid', label='AO')

        axis.legend(handles=[ts, hur_cat1, hur_cat2, hur_cat3,
                        hur_cat4, bt, awo, awo_ws], 
                        prop = { "size": fontsize }, bbox_to_anchor=(xlocator, ylocator),
                        frameon=True, fancybox=True, shadow=True, loc=loc)
        

        
def Track_Legend_v2(axis, fontsize, markersize, xlocator, ylocator, loc):
        #Handles for legend
        bt = mlines.Line2D([], [], color='black', marker='s', markersize=markersize, ls='solid', label='best track')
        awo = mlines.Line2D([], [], color='red', marker='o', markersize=markersize, ls='solid', label='$AWO$-$CTL$')
        awo_ws = mlines.Line2D([], [], color='cyan', marker='*', markersize=markersize,ls='solid', label='$AWO_{ws}$-$EXP$')
        axis.legend(handles=[bt, awo_ws, awo], 
                        prop = { "size": fontsize }, bbox_to_anchor=(xlocator, ylocator),
                        frameon=True, fancybox=True, shadow=True, loc=loc)
        
def add_corner_label(ax, text, fontsize):
    ax.text(0.03, 0.92, text, transform=ax.transAxes, bbox=dict(facecolor='darkgrey', alpha=0.8), fontsize=fontsize, fontweight='bold')



def add_sst_colorbar(axis, data, ticks):
        cbaxes = inset_axes(axis, width="8%", height="80%", loc='upper right')
        cb = plt.colorbar(data, cax=cbaxes, ticks=ticks, orientation='vertical')
        cb.ax.set_title('$^{\circ}C$', color='red', fontsize=10)
        cbaxes.tick_params(labelsize=10, colors='red', bottom=False, top=True, labeltop=True, labelbottom=False)

def add_colorbar(axis, data, ticks, color, units):
        cbaxes = inset_axes(axis, width="8%", height="80%", loc='upper right')
        cb = plt.colorbar(data, cax=cbaxes, ticks=ticks, orientation='vertical')
        cb.ax.set_title(units, color=color, fontsize=10)
        cbaxes.tick_params(labelsize=10, colors=color, bottom=False, top=True, labeltop=True, labelbottom=False)

def add_surf_curr_colorbar(axis, data, ticks):
        cbaxes = inset_axes(axis, width="8%", height="80%", loc='upper right')
        cb = plt.colorbar(data, cax=cbaxes, ticks=ticks, orientation='vertical')
        cb.ax.set_title('$m/s$', color='red', fontsize=10)
        cbaxes.tick_params(labelsize=10, colors='red', bottom=False, top=True, labeltop=True, labelbottom=False)



def add_storm_relative_essentials(axis):
    # axis.axvline(0, color='k', linestyle='solid', linewidth=0.5)
    # axis.axhline(0, color='k', linestyle='solid', linewidth=0.5)
    # axis.arrow(0, 0, 0, 250, width=1.5, color='k')

    axis.axvline(0, color='k', linestyle='solid', linewidth=0.5)
    axis.axhline(0, color='k', linestyle='solid', linewidth=0.5)
    axis.arrow(0, 0, 0, 175, width=0.75, color='k')
    first_circle_wspd = plt.Circle( (0, 0), 50, lw=0.5, fill = False, color='k')
    second_circle_wspd = plt.Circle( (0, 0), 150, lw=0.5, fill = False, color='k')
    third_circle_wspd = plt.Circle( (0, 0), 250, lw=0.5, fill = False, color='k')
    axis.add_artist(first_circle_wspd)
    axis.add_artist(second_circle_wspd)
    axis.add_artist(third_circle_wspd)

def cross_section_essential(x_pos, y_pos, fontsize, lw):
    #Vertical line separating the right and left of the storm
    plt.axvline(x=0, color = 'k', linestyle='--', lw=lw)
    plt.text(x_pos,y_pos, 'R', fontweight='bold', fontsize=fontsize)
    plt.text(-1 * x_pos,y_pos, 'L', fontweight='bold', fontsize=fontsize)

def add_storm_relative_axis_labels(axis, fontsize, labelpad, labelsize, xticks, yticks):
    axis.set_ylabel('Distance to Storm Center (km)', fontsize=fontsize, labelpad=labelpad)
    axis.set_xlabel('Distance to Storm Center (km)', fontsize=fontsize, labelpad=labelpad)
    axis.set_xticks(ticks=xticks)
    axis.set_yticks(ticks=yticks)
    axis.tick_params(axis='both', which='major', labelsize=labelsize, length=1, width=1, pad=0.5)

def add_vertical_cross_section_axis_labels(axis, fontsize, labelpad, labelsize, xticks, yticks):
    axis.set_ylabel('Depth (m)', fontsize=fontsize, labelpad=labelpad)
    axis.set_xlabel('Distance to Storm Center (km)', fontsize=fontsize, labelpad=labelpad)
    axis.set_xticks(ticks=xticks)
    axis.set_yticks(ticks=yticks)
    axis.tick_params(axis='both', which='major', labelsize=labelsize, length=1, width=1, pad=0.5)

def add_corner_label(ax, x_pos, y_pos, text, fontsize):
    ax.text(x_pos, y_pos, text, transform=ax.transAxes, bbox=dict(facecolor='darkgrey', alpha=0.8), fontsize=fontsize, fontweight='bold')

def get_mixed_layer_depth(temp_data, mld_data, Y_vals):
    for ii in range(len(mld_data)):
            for kk in range(11,500): # MLD not likely deeper than 500 m.
                  tdiff = temp_data[kk,ii] - temp_data[10,ii]
                  if tdiff < -0.5:
                        break
            mld_data[ii] = Y_vals[kk]
    return mld_data;
# Old fontsize =9

def extract_data_ahead_of_storm(data1, data2, n):
    data_diff = data1 - data2

    #Empty list to store data after removing nans
    cleaned_data = []

    #Empty list to store all the median values for each time step
    median_of_data = create_nan_list(n) 

    for i in range(n):
        # print(i)
        data_diff_qc = xr.where(data_diff[i].y > 0 , data_diff[i], np.nan)

        #Converting data to list
        data_with_nans = data_diff_qc.values.flatten().tolist()

        #Removing all nans because whisker plots will not work with nans
        data_qc = [value for value in data_with_nans if not math.isnan(value)]

        #Computing and storing the median of each period
        median_of_data[i] = np.nanmedian(data_qc)

        #Combining all list into a single one
        cleaned_data.append(data_qc)

    return cleaned_data, median_of_data;

def create_nan_list(num_nans):
    return [float('nan')] * num_nans

def extract_data_behind_of_storm(data1, data2, n):
    data_diff = data1 - data2

    #Empty list to store data after removing nans
    cleaned_data = []

    #Empty list to store all the median values for each time step
    median_of_data = create_nan_list(n) 

    for i in range(n):
        # print(i)
        data_diff_qc = xr.where(data_diff[i].y < 0 , data_diff[i], np.nan)

        #Converting data to list
        data_with_nans = data_diff_qc.values.flatten().tolist()

        #Removing all nans because whisker plots will not work with nans
        data_qc = [value for value in data_with_nans if not math.isnan(value)]

        #Computing and storing the median of each period
        median_of_data[i] = np.nanmedian(data_qc)

        #Combining all list into a single one
        cleaned_data.append(data_qc)

    return cleaned_data, median_of_data;

def generate_best_fit_line(x_data, y_data, polfit):
    coeff = np.polyfit(x_data, y_data, polfit)
    slope = coeff[0]
    intercept = coeff[1]

    bf_line = slope * x_data + intercept

    return bf_line, slope, intercept;

def extract_masked_data_of_interest(data1, data2, index_of_mask, n):

    data_diff = data1 - data2
    #Empty list to store data after removing nans
    cleaned_data = []

    #Empty list to store all the median values for each time step
    median_of_data = create_nan_list(n) 

    for i in range(8):
        #Masking data
        data_array = data_diff[i].to_masked_array()
        data_array[index_of_mask] = np.nan
        data_masked_qc = data_array.data

        #Converting data to list
        data_with_nans = data_masked_qc.flatten().tolist()

        #Removing all nans because whisker plots will not work with nans
        data_qc = [value for value in data_with_nans if not math.isnan(value)]

        #Computing and storing the median of each period
        median_of_data[i] = np.nanmedian(data_qc)

        #Combining all list into a single one
        cleaned_data.append(data_qc)

        # print(median_of_data)

    return cleaned_data, median_of_data;

def add_axis_label_to_dub_axis_summary_plot(axis, twin_axis, ylabel1, ylabel2, fontsize, labelpad, labelsize, color1, color2, rotation, ymin, ymax, length, width):
    axis.set_ylabel(ylabel1, fontsize=fontsize, labelpad=labelpad, color=color1)
    twin_axis.set_ylabel(ylabel2, fontsize=fontsize, labelpad=labelpad, color=color2)
    axis.set_xlabel('Date', fontsize=fontsize, labelpad=labelpad)
    axis.set_ylim([ymin,ymax])
    twin_axis.set_ylim([ymin,ymax])

    axis.tick_params(axis='both', which='major', labelsize=labelsize)
    axis.tick_params(axis='x', which='major', labelrotation=rotation)
    axis.tick_params(axis='y', which='major', colors=color1)
    twin_axis.tick_params(axis='y', which='major', colors=color2, labelsize=labelsize, length=length, width=width, pad=0.5)
    axis.tick_params(axis='both', which='major', labelsize=labelsize, length=length, width=width, pad=0.5)

def add_axis_label_to_single_axis_summary_plot(axis, ylabel, fontsize, labelpad, labelsize, rotation, ymin, ymax, length, width):
    axis.set_ylabel(ylabel, fontsize=fontsize, labelpad=labelpad)
    axis.set_xlabel('Date', fontsize=fontsize, labelpad=labelpad)
    axis.set_ylim([ymin,ymax])
    axis.tick_params(axis='both', which='major', labelsize=labelsize)
    axis.tick_params(axis='x', which='major', labelrotation=rotation)
    axis.tick_params(axis='both', which='major', labelsize=labelsize, length=length, width=width, pad=0.5)

def add_legend_to_summary_plot(axis, color1, color2, label1, label2, loc, size):
    swh_med= mlines.Line2D([], [], color=color2, ls='solid', label='$trend$')
    wspd_med = mlines.Line2D([], [], color=color1, ls='solid', label='$trend$')
    swh_fit= mlines.Line2D([], [], color=color1, markersize=1, marker='o', ls='', label=label1)
    wspd_fit= mlines.Line2D([], [], color=color2, markersize=1, marker='o', ls='', label=label2)
    best_fit = mlines.Line2D([], [], color='blue',ls='dashed', label='$y=0$')

    axis.legend(handles=[swh_med, wspd_med, wspd_fit, swh_fit], 
            ncol=2, loc=loc, prop = {"size": size}, frameon=True, fancybox=True, shadow=True)


