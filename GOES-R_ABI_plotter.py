
#-------------------------------------------
#
# FILENAME:	GOES-R_ABI_plotter.py
# 		
# CREATED: 	01.04.2022 - dserke
#
#
#-------------------------------------------

#-------------------------------------------
# IMPORT LIBRARIES
#-------------------------------------------
#Import libraries and settings
#Library to perform array operations
import numpy              as     np 

#Libraries for making plots
import matplotlib         as     mpl
from   matplotlib         import pyplot as plt
import matplotlib.ticker  as     ticker

#Libaries for drawing maps
import cartopy
from   cartopy            import crs as ccrs
import cartopy.feature    as     cfeature
from   cartopy.feature    import NaturalEarthFeature
from   cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#Library for creating animations
from   PIL                import Image

#Library for accessing files in the directory
import os

#Library to read in netCDF files
from   netCDF4            import Dataset

#Library for using math functions
import math

#Library for collecting lists of files from folders
import glob

import warnings
warnings.filterwarnings('ignore')

#Sets font size to 12
plt.rcParams.update({'font.size': 12})

#Option to keep numpy from printing in scientific notation by default
np.set_printoptions(suppress = True)

#-------------------------------------------
# DEFINE CONSTANTS
#-------------------------------------------
# input data file dir
ABI_file_path    = '/d1/serke/projects/PYRO_detect_ATEC/data/ABI/nc/'

# input data file name
ABI_file_name    = 'OR_ABI-L2-FDCC-M6_G16_s20213641901172_e20213641903545_c20213641904147.nc'

# output image path
ABI_image_path   = '/d1/serke/projects/PYRO_detect_ATEC/images/20211230_KFTG/ABI/'

# Define flags
info_to_GUI_flag = False

#-------------------------------------------
# DEFINE FUNCTIONS
#-------------------------------------------
# ALGORITHM-RELATED FUNCTIONS
# ... Algorithm to convert latitude and longitude radian values to degrees (as a callable function)
def Degrees(file_id):
    proj_info        = file_id.variables['goes_imager_projection']
    lon_origin       = proj_info.longitude_of_projection_origin
    H                = proj_info.perspective_point_height+proj_info.semi_major_axis
    r_eq             = proj_info.semi_major_axis
    r_pol            = proj_info.semi_minor_axis 
    #Data info
    lat_rad_1d       = file_id.variables['x'][:]
    lon_rad_1d       = file_id.variables['y'][:]
    #Create meshgrid filled with radian angles
    lat_rad, lon_rad = np.meshgrid(lat_rad_1d,lon_rad_1d)
    #lat/lon calculus routine from satellite radian angle vectors
    lambda_0         = (lon_origin*np.pi)/180.0
    a_var            = np.power(np.sin(lat_rad),2.0) + (np.power(np.cos(lat_rad),2.0)*(np.power(np.cos(lon_rad),2.0)
    +(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(lon_rad),2.0))))
    b_var            = -2.0*H*np.cos(lat_rad)*np.cos(lon_rad)
    c_var            = (H**2.0)-(r_eq**2.0)
    r_s              = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)
    s_x              = r_s*np.cos(lat_rad)*np.cos(lon_rad)
    s_y              = - r_s*np.sin(lat_rad)
    s_z              = r_s*np.cos(lat_rad)*np.sin(lon_rad)
    Lat              = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
    Lon              = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)
    return Lat, Lon

# ... Select and process Mask data from a single file (as a callable function)
def Mask_Data(file_id):
    #Read in Power data
    Mask_data        = file_id.variables['Mask'][:,:]
    #Select quality of AOD data pixels using the "DQF" variable
    #High quality: DQF = 0, Medium quality: DQF = 1, Low quality: DQF = 2, not retrieved (NR): DQF = 3
    #Science team recommends using High and Medium qualities for operational applications- 
    #(e.g.,mask low quality and NR pixels)
    DQF              = file_id.variables['DQF'][:,:]
    Quality_Mask     = (DQF > 1)
    Mask             = np.ma.masked_where(Quality_Mask, Mask_data)
    return Mask

# ... Select and process Power data from a single file (as a callable function)
def Power_Data(file_id):
    #Read in Power data
    Power_data       = file_id.variables['Power'][:,:]
    #Select quality of AOD data pixels using the "DQF" variable
    #High quality: DQF = 0, Medium quality: DQF = 1, Low quality: DQF = 2, not retrieved (NR): DQF = 3
    #Science team recommends using High and Medium qualities for operational applications- 
    #(e.g.,mask low quality and NR pixels)
    DQF              = file_id.variables['DQF'][:,:]
    Quality_Mask     = (DQF > 1)
    Power             = np.ma.masked_where(Quality_Mask, Power_data)
    return Power

# ... Select and process DQF data from a single file (as a callable function)
def DQF_Data(file_id):
    #Read in AOD data
    DQF              = file_id.variables['DQF'][:,:]
    return DQF

# PLOTTING-RELATED FUNCTIONS
# ... Plotting settings for Mask data (as a callable function)
def Mask_Data_Settings():
    #Create custom continuous colormap for Mask data
    #.set_over sets color for plotting data > max
    color_map = mpl.colors.LinearSegmentedColormap.from_list('custom_Mask', [(0, 'indigo'),(.10, 'mediumblue'), 
                                                             (.20, 'blue'), (.30, 'royalblue'), (.40, 'skyblue'), 
                                                             (.50, 'cyan'), (.60, 'yellow'), (.70, 'orange'), 
                                                             (.80, 'darkorange'), (.90, 'red'), (1.00, 'firebrick')], N = 150)
    color_map.set_over('darkred')
    #Set range for plotting Mask data (data min, data max, contour interval) (MODIFY contour interval)
    #interval: 0.1 = runs faster/coarser resolution, 0.01 = runs slower/higher resolution
    data_range = np.arange(0, 11, 0.05)
    return color_map, data_range
  
# ... Create Mask colorbar, independent of plotted data (as a callable function)
# ...... cbar_ax are dummy variables
# ...... Location/dimensions of colorbar set by .set_position (x0, y0, width, height) to scale automatically with plot
def Mask_Colorbar():
    last_axes = plt.gca()
    cbar_ax   = fig.add_axes([0, 0, 0, 0])
    plt.draw()
    posn      = ax.get_position()
    cbar_ax.set_position([0.35, posn.y0 - 0.07, 0.3, 0.02])
    color_map = mpl.colors.LinearSegmentedColormap.from_list('custom_Mask', [(0, 'indigo'),(.10, 'mediumblue'), 
                                                             (.20, 'blue'), (.30, 'royalblue'), (.40, 'skyblue'), (.50, 'cyan'), (.60, 'yellow'), (.70, 'orange'), 
                                                             (.80, 'darkorange'), (.90, 'red'), (1.00, 'firebrick')], N = 150)
    color_map.set_over('darkred')
    norm      = mpl.colors.Normalize(vmin = 0, vmax = 1.00)
    cb        = mpl.colorbar.ColorbarBase(cbar_ax, cmap = color_map, norm = norm, orientation = 'horizontal', 
    ticks     = [0, .25, .50, .75, 1.00], extend = 'max')
    cb.set_label(label = 'Mask', size = 'medium', weight = 'bold')
    cb.ax.set_xticklabels(['0', '.25', '.50', '.75', '1.00'])
    cb.ax.tick_params(labelsize = 'medium')
    plt.sca(last_axes)
    
# ... Format map with Plate Carree projection (as a callable function)
def ABI_Map_Settings_PC(ax):
    #Set up and label the lat/lon grid
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_xticks([-160, -140, -120, -100, -80, -60, -40, -20], crs = ccrs.PlateCarree())
    ax.set_yticks([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80], crs = ccrs.PlateCarree())  
    #Set lat/lon ticks and gridlines
    ax.tick_params(length = 0)
    ax.grid(linewidth = 0.5, zorder = 3)
    #Draw coastlines/borders using Cartopy; zorder sets drawing order for layers
    ax.coastlines(resolution = '50m', zorder = 3)
    ax.add_feature(cfeature.BORDERS, zorder = 3)
    ax.add_feature(cfeature.NaturalEarthFeature(category = 'cultural', name = 'admin_1_states_provinces', 
                                                scale = '50m'), facecolor = 'none', lw = 0.5, edgecolor = 'black', zorder = 4)
    ax.add_feature(cfeature.NaturalEarthFeature(category = 'cultural', name = 'urban_areas', 
                                                scale = '50m'), facecolor = 'none', lw = 0.8, edgecolor = 'red', zorder = 4)
    ax.add_feature(cfeature.NaturalEarthFeature(category = 'physical', name = 'ocean', scale = '50m'), 
                                                facecolor = 'lightgrey')
    ax.add_feature(cfeature.NaturalEarthFeature(category = 'physical', name = 'land', scale = '50m'), 
                                                facecolor = 'grey')
    ax.add_feature(cfeature.NaturalEarthFeature(category = 'physical', name = 'lakes', scale = '50m'), 
                                                facecolor = 'lightgrey', edgecolor = 'black', zorder = 2)
    #Set domain for map [x0, x1, y0, y1] 
    #  Default longitude extent (x0, x1) values: G16 = (-135, -65); G17 = (-170, -100)
    #  Use 180 degrees for longitude coordinates (i.e, -100 = 100 degrees W)
    #  NOTE: comment out (add leading ##) the line below to automatically set domain to extend to limits of data
    ax.set_extent([-105.7, -103.0, 39.4, 42.0], crs       = ccrs.PlateCarree())
    ##ax.set_extent([-135, -65, 15, 55], crs       = ccrs.PlateCarree())

#-------------------------------------------
# DEFINE INPUT FILES
#-------------------------------------------
# Define single file name
fname            = ABI_file_path+ABI_file_name
# Define file lists
file_list        = sorted(glob.glob(ABI_file_path+'/*FDCC*.nc'))

#Loop through data files, making/saving a figure for each data file
for x in file_list:

    #-------------------------------------------
    # LOAD INPUT FIELDS FROM INPUT FILE
    #-------------------------------------------
    # Set the file name to read
    file_id          = Dataset(x)
    # Select and process FDCC data fields
    Mask             = Mask_Data(file_id)
    Power            = Power_Data(file_id)
    DQF              = DQF_Data(file_id)
    # Read in latitutude and longitude values in degrees
    Lat, Lon         = Degrees(file_id)

    #-------------------------------------------
    # DISPLAY/GET INFO FROM INPUT DATA
    #-------------------------------------------
    if info_to_GUI_flag == True:
        #Check the contents of the entire file
        print(file_id)
        #Check the AOD variable metadata
        print(file_id.variables['Mask'])
        print(file_id.variables['Power'])
        print(file_id.variables['DQF'])
        #Check the AOD array values
        print(file_id.variables['Mask'][:,:])
        print(file_id.variables['Power'][:,:])
        print(file_id.variables['DQF'][:,:])
        #Check the spatial resolution of the data
        print((file_id.__getattr__('title')),'spatial resolution is', (file_id.__getattr__('spatial_resolution')))
        #Check the units for the variables of interest (note: "1" means unitless)
        print('Mask unit is', (file_id.variables['Mask'].__getattr__('units')))
        print('Power unit is', (file_id.variables['Power'].__getattr__('units')))
        print('DQF unit is', (file_id.variables['DQF'].__getattr__('units')))
        print('Latitude unit is', (file_id.variables['x'].__getattr__('units')))
        print('Longitude unit is', (file_id.variables['y'].__getattr__('units')))
        #Check the data types for the variables of interest
        print('Mask data type is', (file_id.variables['Mask'][:,:].dtype))
        print('Power data type is', (file_id.variables['Power'][:,:].dtype))
        print('DQF data type is', (file_id.variables['DQF'][:,:].dtype))
        print('Latitude data type is', (file_id.variables['x'][:].dtype))
        print('Longitude data type is', (file_id.variables['y'][:].dtype))
        print('Mask: minimum value is ' + str(np.min(Mask)) + ';' + ' maximum value is ' + str(np.max(Mask)))
        print('Power: minimum value is ' + str(np.min(Power)) + ';' + ' maximum value is ' + str(np.max(Power)))
        print('DQF: minimum value is ' + str(np.min(DQF)) + ';' + ' maximum value is ' + str(np.max(DQF)))
        print('Latitude: minimum value is ' + str(np.min(Lat)) + ' degrees;' + ' maximum value is ' 
              + str(np.max(Lat)) + ' degrees')
        print('Longitude: minimum value is ' + str(np.min(Lon)) + ' degrees;' + ' maximum value is ' 
              + str(np.max(Lon)) + ' degrees')

    #-------------------------------------------
    # PLOTTING
    #-------------------------------------------
    # ... Plotting settings for data
    color_map, data_range = Mask_Data_Settings()

    # PLOT 'Mask' DATA - CONUS VIEW
    # ... default central_longitude = 0
    fig                   = plt.figure(figsize=(24, 30))
    # ... Set up figure and map projection: PlateCarree(central_longitude)
    # ... Plate Carree: equidistant cylindrical projection w/equator as the standard parallel 
    ax                    = fig.add_subplot(1,1,1, projection = ccrs.PlateCarree())
    # ... Format map with Plate Carree projection
    ABI_Map_Settings_PC(ax)
    # ... Add and format title
    # ......Reverse indexing (from right to left) of file name automatically adds satellite, time, and year to title
    plt.title('GOES-' + fname[-53:-51] + '/ABI\nHigh + Medium Quality Fire Mask\n' + x[-42:-40] + ':' + x[-40:-38]
              + ' UTC, 30 May ' + x[-49:-45], y = 1.025, ma = 'center', size = 15, weight = 'bold')
    # ... Add Mask colorbar
    Mask_Colorbar()
    if Mask.count() > 0:
        #Create filled contour plot of data
        Plot      = ax.contourf(Lon, Lat, Mask/100, data_range, cmap = color_map, extend = 'both', zorder = 3, 
        transform = ccrs.PlateCarree())
    else:
        pass
    # ... Show figure
    plt.show()
    # ... Save figure as a .png file
    #  dpi sets the resolution of the digital image in dots per inch
    Mask_img_fname = ABI_image_path+'G16_CONUS_ABI_FIREMASK_'+x[-49:-45]+'MMDD_'+x[-42:-38]+'Z.png'
    fig.savefig(Mask_img_fname, bbox_inches = 'tight', dpi = 150)
    # ... Erase plot so we can build the next one
    plt.close()
    

    # PLOT 'Power' DATA - CONUS VIEW
    # ... default central_longitude = 0
    fig                   = plt.figure(figsize=(24, 30))
    # ... Set up figure and map projection: PlateCarree(central_longitude)
    # ... Plate Carree: equidistant cylindrical projection w/equator as the standard parallel 
    ax                    = fig.add_subplot(1,1,1, projection = ccrs.PlateCarree())
    # ... Format map with Plate Carree projection
    ABI_Map_Settings_PC(ax)
    # ... Add and format title
    # ......Reverse indexing (from right to left) of file name automatically adds satellite, time, and year to title
    plt.title('GOES-' + fname[-53:-51] + '/ABI\nHigh + Medium Quality Fire Power\n' + x[-42:-40] + ':' + x[-40:-38]
              + ' UTC, 30 May ' + x[-49:-45], y = 1.025, ma = 'center', size = 15, weight = 'bold')
    # ... Add Mask colorbar
    Mask_Colorbar()
    # ... Plotting settings for data
    color_map, data_range = Mask_Data_Settings()
    if Power.count() > 0:
        #Create filled contour plot of data
        Plot      = ax.contourf(Lon, Lat, Power/1000, data_range, cmap = color_map, extend = 'both', zorder = 3, 
        transform = ccrs.PlateCarree())
    else:
        pass
    # ... Show figure
    plt.show()
    # ... Save figure as a .png file
    #  dpi sets the resolution of the digital image in dots per inch
    Power_img_fname = ABI_image_path+'G16_CONUS_ABI_FIREPOWER_'+x[-49:-45]+'MMDD_'+x[-42:-38]+'Z.png'
    fig.savefig(Power_img_fname, bbox_inches = 'tight', dpi = 150)
    # ... Erase plot so we can build the next one
    plt.close()
    
# Make an animation of AOD figures using python image library (Pillow)
# ... Pillow is preferred for AOD animations because it retains the features of continuous colorbars relatively well
# ... Collect all of the AOD graphics files (figures) in given subdirectory
Power_img_file_list = sorted(glob.glob(ABI_image_path+'G16_CONUS_ABI_FIREPOWER_'+x[-49:-45]+'MMDD_*Z.png'))
Mask_img_file_list  = sorted(glob.glob(ABI_image_path+'G16_CONUS_ABI_FIREMASK_'+x[-49:-45]+'MMDD_*Z.png'))
# ... Create an empty list to store figures
Power_frames        = []
Mask_frames         = []
# ... Loop through graphics files and append
for x in Power_img_file_list:
    Power_new_frame = Image.open(x)
    Power_frames.append(Power_new_frame)
for x in Mask_img_file_list:
    Mask_new_frame = Image.open(x)
    Mask_frames.append(Mask_new_frame)    
# ... Save animation
#Find the saved animation in your "current working directory" (folder containing Python Notebook file,- 
#netCDF-4 data files)
#Duration is speed of frame animation in ms (e.g., 1000 ms = 1 second between frames)
#Loop sets time before animation restarts (e.g., loop = 0 means animation loops continuously with no delay)
Power_frames[0].save('G16_CONUS_ABI_FIREPOWER-Animation.gif', format = 'GIF', append_images = Power_frames[1:], save_all = True, duration = 1000, loop = 0)
Mask_frames[0].save('G16_CONUS_ABI_FIREMASK-Animation.gif',   format = 'GIF', append_images = Mask_frames[1:],  save_all = True, duration = 1000, loop = 0)
# ... Close the graphics files we opened
for x in Power_img_file_list:
    Power_new_frame.close()
for x in Mask_img_file_list:
    Power_new_frame.close()    
print('Animation done!')