#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 17 2019
Code to open CWATM diagnostics for individual monthly or daily files 
and create lat-lon timeseries for each variable and save as netcdf.
@author: M. Puma, N. Woodruff
"""

# Import Modules and define functions
import numpy as np
import netCDF4

# For plotting a rectangle on the maps
def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys,latlon = True, color='k', linestyle='--', linewidth=3)

## Directory Names 
dir_name = '/Users/puma/Soil_TerraE_Diag/'            # for Inputs
out_dir = '/Users/puma/Soil_TerraE_Diag/processed/'   # for Outputs

### Load time invariant data, needed for some calculations.
# Pull out time invariant variables (gridcell area, lat, lon)
fname = '/Users/puma/Soil_TerraE_Diag/fixed/fixedvar.modelE.nc'
            
# Month Vectors
mons     = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
mons_txt = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

# Set plot styles
# Formatting for titles
fontdict_title = {'fontsize': 36}
fig_size = np.array([10,10])

# Formatting for figures
style_new = {'xtick.direction': 'in', \
             'ytick.direction': 'in', \
             'font.sans-serif': 'Arial'}

# Open this netcdf file
ncfile = netCDF4.Dataset(fname) # Load Dimension/Fixed Variables

# Load variables
lat      = ncfile.variables['lat'][:]         # latitude, degrees
lon      = ncfile.variables['lon'][:]         # longitude, degrees

# Close File
ncfile.close

### Setup for Variables/Simulation/Etc
# Choose simulation
sim_name='P01CarbonSHPs'; scen='CarbonSHPs';      memb='P01'; memb_new='D01';
numstr='01'; sim_name='P'+numstr+'CarbonSHPs'; scen='CarbonSHPs'; memb='P'+numstr; memb_new='D'+numstr; 
yrs_split = [np.arange(1860,1869+1)]

# Split up Years Into Different Blocks ot time
# FOR HISTORICAL RUNS
#yrs_split = [  np.arange(1870,2014+1)]
#yrs_split = [  np.arange(2015,2114+1)]

# FOR PI CTRL RUNS
#sim_name='E01PICFO';    scen='piCtrlFixSST'; memb='E01'; memb_new='D01';
#yrs_split = [  np.arange(1000,1500+1) ]

#sim_name='E01PICQF';    scen='piCtrlQFLX'; memb='E01'; memb_new='D01';
#yrs_split = [  np.arange(1000,1500+1) ]

# Variable and output information
# variable name, new name, input file type, units, long name, other info (atmo vs land, monthly vs daily, etc)
var_meta = [   ['Tatm',  'Tatm',   'aij',  'C',     'atmospheric temperature',  'Amon'], \
               ['tsurf', 'ta',     'aij',  'C',     'surface air temperature',  'Amon'], \
               ['prec',  'pr',     'aij',  'mm/day', 'precipitation',           'Amon'], \
               ['sst',   'sst',    'aij',  'C',     'sea surface temperature',  'Amon'], \
               ['evap',  'evap',   'aij',  'mm/day','evaporation',              'Lmon'], \
               ['irrig_w','irrig_w','aij', 'mm/day','total irrigation',         'Lmon'], \
               ['irrig_gw','irrig_gw','aij', 'mm/day','fossil irrigation',      'Lmon'], \
            ] 

### Concatenate output over time and save variables into separate netcdf files

# Time vector (months since )
time_name  = 'time'

# Loops: Variable, Year, Month----------------------------------------------------------------------------------------

# Number of variables
num_var = np.shape(var_meta)[0]
num_blk = np.shape(yrs_split)[0]

# Time interval loop (do TEN YEAR BLOCKS)
for i_time in np.arange(num_blk):
    
    # Pull out time information
    yrs        = yrs_split[i_time][:]
    time_units = 'months since '+np.str(yrs[0])+'-1-1 0:0:0'
    time_vect  = np.arange(0.5,yrs.size*12)
    
    #-----------------------------------------------------------------------------------------------------------------
        
    # Variable Loop
    for i_var in np.arange(num_var):
    
        # Pull out variable information
        var_name   = var_meta[i_var][0]
        cf_var     = var_meta[i_var][1] 
        file_form  = var_meta[i_var][2] 
        units      = var_meta[i_var][3]
        long_name  = var_meta[i_var][4]
        other_info = var_meta[i_var][5]
    
        # Initialize storage array (first slice is NaN; will remove later)
        var_all_time = np.zeros((1,lat.size,lon.size))*np.nan
    
        # Year Loop
        for i_yr in enumerate(yrs):
    
            # Month Loop    
            for i_mon in enumerate(mons_txt):
        
                # Create file name for current file
                fname = dir_name+i_mon[1]+np.str(i_yr[1])+'.'+file_form+sim_name+'.nc'
                
                # Open this netcdf file
                ncfile = netCDF4.Dataset(fname)   

                # Initialize storage array and loop through each layer
                store_var = np.zeros((1,lat.size,lon.size))*np.nan
                                    
                # Load data
                var  = np.array(ncfile.variables[var_name][:]);      

                # Store in array 
                store_var[0,:,:] = var                   
                    
                # Concatenate in time
                var_all_time = np.concatenate((var_all_time,store_var),axis=0)
                
                # Close NetCDF File
                ncfile.close()
                
        #-------------------------------------------------------------------------------------------------------------
        # Trim off placeholder at beginning
        var_all_time = var_all_time[1:,:,:]
        
        # Save to netcdf file
        
        # Output filename and number of layers
        out_file  = cf_var+'_'+other_info+'_ModelE'+'_'+scen+'_'+memb_new+'_'+np.str(np.min(yrs))+'01'+np.str(np.max(yrs))+'12.nc'
        fname_out = out_dir+out_file
        print(var_name+': '+fname_out)
        
        # Create and open output file
        ncout = netCDF4.Dataset(fname_out,'w',clobber=True,format="NETCDF4_CLASSIC")

        # Create Dimensions
        ncout.createDimension('lat',   np.size(lat))
        ncout.createDimension('lon',   np.size(lon))
        ncout.createDimension('time',  np.size(time_vect))

        # Create Variables (info: name, precision, dimensions)
        # dimensions
        lat_nc   = ncout.createVariable('lat', float, ('lat'), zlib=True)
        lon_nc   = ncout.createVariable('lon', float, ('lon'), zlib=True)
        time_nc  = ncout.createVariable('time', float, ('time'), zlib=True)
    
        # Two-Dimensional Variables
        var_nc  = ncout.createVariable(cf_var, float, ('time','lat','lon'), zlib=True)

        # Write out dimensional data
        lat_nc[:]  = lat;         lat_nc.long_name    = 'latitude';       lat_nc.units   = 'degrees_north';
        lon_nc[:]  = lon;         lon_nc.long_name    = 'longitude';      lon_nc.units   = 'degrees_east';
        time_nc[:] = time_vect;   time_nc.long_name   = time_name;        time_nc.units  = time_units;   

        # Write out data
        var_nc[:] = var_all_time; var_nc.long_name = long_name;  var_nc.units = units;

        # File Information
        ncout.comment = 'Variable from CWATM experiment ('+sim_name+')'

        # Close the file
        ncout.close()