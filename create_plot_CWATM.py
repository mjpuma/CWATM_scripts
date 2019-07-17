#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 17 2019
Code to analyze and plot CWATM diagnostics processed by script 
"create_timeseries_CWATM.py"
@author:  M. Puma; N. Woodruff
"""

# Import Modules and define functions
import calendar
import datetime
import os
import numpy as np
import numpy.ma as ma
import netCDF4
import matplotlib
import copy
from matplotlib import pyplot as plt
import scipy
import scipy.signal
import scipy.io as sio
import seaborn as sns
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from IPython.display import display
#from mpl_toolkits.basemap import Basemap, cm, maskoceans

# cartopy stuff
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry.polygon import LinearRing
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Borders for mapping: Cultural borders
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
# Coastline
newcoast = cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='k',
                                        facecolor='none')
#Lakes
newlake = cfeature.NaturalEarthFeature('physical', 'lakes', '10m',
                                        edgecolor='k',
                                        facecolor='none')


# For plotting a rectangle on the maps
#def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax):
#    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
#    ys = [latmin,latmin,latmax,latmax,latmin]
#    bmap.plot(xs, ys,latlon = True, color='k', linestyle='--', linewidth=3)
def plot_rectangle(ax, lonmin,lonmax,latmin,latmax):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    #ax.plot(xs, ys,latlon = True, color='k', linestyle='--', linewidth=3)
    ax.plot(xs,ys,color='k',linestyle='--',linewidth=3,transform=ccrs.PlateCarree())


## Directories
dir_name  = '/Users/puma/Soil_TerraE_Diag/processed';

# Month Vector
mons     = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

# Set plot styles
# Formatting for titles
fontdict_title = {'fontsize': 36}
fig_size = np.array([10,10])

# Formatting for figures
style_new = {'xtick.direction': 'in', \
             'ytick.direction': 'in', \
             'font.sans-serif': 'Arial'}

# Years to Analyze (including calculation of climatology baseline)
yrs_clim = np.arange(1900,1929+1)
yrs_hist = np.arange(1948,1957+1)
yrs_fut  = np.arange(2048,2057+1)

# Regions for averaging (southwest & southeast)
latlon_swusa = np.array([-115,-100,28,40])
latlon_seusa = np.array([-95,-75,28,37])

# Indices for different seasonal intervals
i_JJA    = np.where((mons>=6) & (mons<=8))[0]
i_OND    = np.where((mons>=10) & (mons<=12))[0]
i_JFM    = np.where((mons>=1) & (mons<=3))[0]
i_AMJJAS = np.where((mons>=4) & (mons<=9))[0]

# Pull out time invariant variables (gridcell area, lat, lon)
fname = '/Users/puma/Soil_TerraE_Diag/fixed/fixedvar.modelE.nc'
            
# Open this netcdf file
ncfile = netCDF4.Dataset(fname) # Load Dimension/Fixed Variables

# Load variables
lat      = ncfile.variables['lat'][:]         # latitude, degrees
lon      = ncfile.variables['lon'][:]         # longitude, degrees
areacell = ncfile.variables['axyp'][:]        # gridcell area, m2
fracland = ncfile.variables['frac_land'][:]   # fraction of land in gridcell (%)
arealand = areacell*(fracland/100)

# Close the netcdf file
ncfile.close 

plt.pcolor(lon,lat,fracland,cmap='viridis')
plt.colorbar()

# Lat/lon Indices for different regions
i_lat_nino12 = np.where( (lat<=0) & (lat>=-10) )[0] ; i_lon_nino12 = np.where( (lon<=-80)  & (lon>=-90) )[0]
i_lat_nino34 = np.where( (lat<=5) & (lat>=-5) )[0] ;  i_lon_nino34 = np.where( (lon<=-120) & (lon>=-170) )[0]
i_lat_nino4A = np.where( (lat<=5) & (lat>=-5) )[0] ;  i_lon_nino4A = np.where( (lon<=180)  & (lon>=160) )[0]
i_lat_nino4B = np.where( (lat<=5) & (lat>=-5) )[0] ;  i_lon_nino4B = np.where( (lon<=-150) & (lon>=-180) )[0]

## Some simple mapping to familiarize myself with cartopy
#* Default boundaries:
#* #ax.add_feature(cfeature.LAKES, linewidth=1, linestyle='-', edgecolor='k')
#* #ax.add_feature(cfeature.COASTLINE, linewidth=1.5, linestyle='-', edgecolor='k')
#* #ax.add_feature(cfeature.STATES, linewidth=.5)
#* #ax.add_feature(countries, linewidth=0.5, linestyle='-', edgecolor='k')
# Setup some map stuff
clevs       = np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
ocean_color = np.float64([209,230,241])/255

# Lat/Lon range for regional map
extent_lonlat = (-135, -60, 15, 65)

# Create New Colormap With a Limited Number of entries
nmap=plt.cm.get_cmap(name='viridis',lut=9) # only needed to set bins for pcolor/pcolormesh
#nmap=plt.cm.get_cmap(name=plt.cm.BrBG,lut=9) # only needed to set bins for pcolor/pcolormesh
#nmap=plt.cm.get_cmap(name='BrBG',lut=12)

# Set Data
#data_map = copy.deepcopy(arealand); data_map = (data_map/1e6) # convert to km2
data_map = fracland/100;

# Mask Ocean Areas
data_map[np.where(fracland<=0)]=np.nan

#---------------------------------------------------------------------------------------------------------------------

# Global Map

# Create Figure and set projection
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
#ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()

#---------------------------------------------------------------------------------------------------------------------

# Create Figure and set projection
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# Plot Data
m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
# Lat/Lon Axis Labels
ax.set_xticks(np.arange(-180,190,10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(-90,100,10), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labelsize=18,direction="in")
ax.set_global()
# Change geographic extent and add cultural/physical features
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=1, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
# Add a Rectangle and some annotation
#plot_rectangle(ax,-120,-90,20,40)
plot_rectangle(ax,-115,-100,28,40)
plot_rectangle(ax,-95,-75,28,37)
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.06,pad=0.06)
cbar.ax.tick_params(labelsize=24)
plt.show()
    

#---------------------------------------------------------------------------------------------------------------------
extent_lonlat = (-135, -65, 15, 65)

# Regional Map, different projection
# Create Figure and set projection
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
#m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
#              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0)
# Add a Rectangle and some annotation
#plot_rectangle(ax,-120,-90,20,40)
plot_rectangle(ax,-115,-100,28,40)
plot_rectangle(ax,-95,-75,28,37)
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()

#---------------------------------------------------------------------------------------------------------------------
extent_lonlat = (-135, -65, 15, 65)

# Regional Map, different projection
# Create Figure and set projection
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0)
# Add a Rectangle and some annotation
#plot_rectangle(ax,-120,-90,20,40)
plot_rectangle(ax,-115,-100,28,40)
plot_rectangle(ax,-95,-75,28,37)
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()

## Calculate SST Climatologies and map anomalies for 1860-1864 and 1865-1869


## Load and calculate global annual temperature, SSTs, and irrigation from historical simulations 
# List of ensemble members
ens_list = ['D01']
#ens_list = ['D01','D02','D03','D04','D05','D06','D07','D08','D09','D10' \
#           ,'D11','D12','D13','D14','D15','D16','D17','D18','D19','D20']

# Unique yrs
yrs = np.arange(1860,1869+1)

# Full year/month vectors
yrs_vect  = np.repeat(np.arange(1860,1869+1),12,axis=0)
mons_vect = np.ndarray.flatten(np.repeat([np.arange(1,12+1)],np.unique(yrs_vect).size,axis=0))

# Storage arrays for time series
tsurf_globe_histfut = np.zeros((np.size(ens_list),yrs.size,12))*np.nan
nino12_histfut = np.zeros((np.size(ens_list),yrs.size,12))*np.nan
nino34_histfut = np.zeros((np.size(ens_list),yrs.size,12))*np.nan
nino4_histfut  = np.zeros((np.size(ens_list),yrs.size,12))*np.nan


# Loop Through Each Ensemble------------------------------------------------------------------------------------------
for i_ens in enumerate(ens_list):
       
    # Load information from current ensemble member
    ens_name = ens_list[i_ens[0]]   
    print(ens_name)
    
    # Surface Air Temperature-----------------------------------------------------------------------------------------
    ncfile     = netCDF4.Dataset(dir_name+'/ta_Amon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')   
    tsurf_all = ncfile.variables['ta'][:]
    ncfile.close
    
    # SSTs------------------------------------------------------------------------------------------------------------
    ncfile    = netCDF4.Dataset(dir_name+'/sst_Amon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')   
    sst_all  = ncfile.variables['sst'][:]
    ncfile.close
    
    # Nested Year/Month Loops For Spatial Averaging-------------------------------------------------------------------
    for i_yr in enumerate(yrs):
        for i_mon in enumerate(mons):

            # Pull out current year/month of data
            i_loc              = np.where((yrs_vect==i_yr[1]) & (mons_vect==i_mon[1]))[0]
            
            # Global Air temperature Calculation
            curr_tsurf         = np.squeeze(tsurf_all[i_loc,:,:])
            tsurf_globe_histfut[i_ens[0],i_yr[0],i_mon[0]] = np.sum(curr_tsurf*areacell)/np.sum(areacell)

            # Global SSTs, for ENSO Index Calculations
            curr_sst    = np.squeeze(sst_all[i_loc,:,:])
            
            # NINO 3.4 Calculation
            sst_enso    = curr_sst[i_lat_nino34,:][:,i_lon_nino34]
            area_enso   = areacell[i_lat_nino34,:][:,i_lon_nino34]
            nino34_histfut[i_ens[0],i_yr[0],i_mon[0]] = np.sum(sst_enso*area_enso)/np.sum(area_enso)

            # NINO 1+2 Calculation
            sst_enso    = curr_sst[i_lat_nino12,:][:,i_lon_nino12]
            area_enso   = areacell[i_lat_nino12,:][:,i_lon_nino12]
            nino12_histfut[i_ens[0],i_yr[0],i_mon[0]] = np.sum(sst_enso*area_enso)/np.sum(area_enso)

            # NINO 4 Calculation (stitch together area from both sides of dateline)
            sst_ensoA    = curr_sst[i_lat_nino4A,:][:,i_lon_nino4A]
            area_ensoA   = areacell[i_lat_nino4A,:][:,i_lon_nino4A]
            sst_ensoB    = curr_sst[i_lat_nino4B,:][:,i_lon_nino4B]
            area_ensoB   = areacell[i_lat_nino4B,:][:,i_lon_nino4B]
        
            sst_enso4  = np.hstack((sst_ensoA,sst_ensoB))
            area_enso4 = np.hstack((area_ensoA,area_ensoB))
        
            nino4_histfut[i_ens[0],i_yr[0],i_mon[0]] = np.sum(sst_enso4*area_enso4)/np.sum(area_enso4)   

## Save monthly NINO 3.4 SSTs            
## if more than 1 ensemble
##nino34_all = nino34_histfut[1,:,:];
# one ensemble only
nino34_all = nino34_histfut[0,:,:];

df_nino34  = pd.DataFrame(index=yrs,data=nino34_all)

df_nino34.to_csv('/Users/puma/Soil_TerraE_Diag/processed/data/nino34_monthly.csv')

## Calculate SST Climatologies and map anomalies for two periods (referred to 
#    as 20th and 21th century respectively)
# Month Loop
yrs_clim_20th = np.arange(1860,1864+1)
yrs_clim_21st = np.arange(1865,1869+1)

# Storage Arrays for SSTs
clim_sst20th = np.zeros((12,lat.size,lon.size))
clim_sst21st = np.zeros((12,lat.size,lon.size))

#------------------------------------------------------------------------------------------------------------------
# Calculate monthly climatology for baseline period
for i_mon in enumerate(mons):
    
    # Climatology (20th Century)
    i_clim = np.where(  (yrs_vect>=np.min(yrs_clim_20th)) & (yrs_vect<=np.max(yrs_clim_20th)) & \
                          (mons_vect==i_mon[1])  )[0]
        
    # Means
    clim_sst20th[i_mon[0],:,:]        = np.mean(sst_all[i_clim,:,:],axis=0)
 
    # Climatology (21st Century)
    i_clim = np.where(  (yrs_vect>=np.min(yrs_clim_21st)) & (yrs_vect<=np.max(yrs_clim_21st)) & \
                          (mons_vect==i_mon[1])  )[0]
        
    # Means
    clim_sst21st[i_mon[0],:,:]        = np.mean(sst_all[i_clim,:,:],axis=0)
    
#------------------------------------------------------------------------------------------------------------------

# Storage Arrays for SST Anomalies
anom_sst20th = np.zeros((yrs.size,12,lat.size,lon.size))
anom_sst21st = np.zeros((yrs.size,12,lat.size,lon.size))

# Year loop-calculate SST Anomalies
for n_yr in enumerate(yrs):
        
    # Indices (all months, current year)
    i_yr = np.where(yrs_vect==n_yr[1])[0]
    
    # Pull out current year temp/prec/etc
    curryr_sst  = sst_all[i_yr,:,:]

    # SST Anomaly
    anom_sst20th[n_yr[0],:,:,:] = curryr_sst-clim_sst20th
    anom_sst21st[n_yr[0],:,:,:] = curryr_sst-clim_sst21st

plt.pcolormesh(clim_sst21st[1,:,],vmin=0,vmax=34,cmap='viridis'),plt.colorbar()

### SST Trend: linear, DJF average, 2015-2100

# Calculate DJF average SST
yr_djf = yrs[1:]
dec_sst_anom = anom_sst20th[0:-1,11,:,:]
jan_sst_anom = anom_sst20th[1:,0,:,:]
feb_sst_anom = anom_sst20th[1:,1,:,:]
djf_sst_anom = (dec_sst_anom+jan_sst_anom+feb_sst_anom)/3

# Year indices
i_yrs = np.where((yr_djf>=1860) & (yr_djf<=1869))

# Storage array for trend/slope of regression
djf_sst_trend = np.zeros((lat.size,lon.size))*np.nan

# Loop through each grid cell
for i_lat in enumerate(lat):
    for i_lon in enumerate(lon):
        djf_sst_trend[i_lat[0],i_lon[0]] = stats.linregress(yr_djf[i_yrs], \
                djf_sst_anom[i_yrs,i_lat[0],i_lon[0]]).slope

# Setup some map stuff
#clevs       = np.array((-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,.1,.2,.3,.4,.5,.6))/20
#clevs       = np.array((-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,.1,.2,.3,.4,.5,.6))/20
#clevs       = np.array((0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09))
clevs       = np.array((0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06))
#clevs       = np.array((-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6))/60
#clevs       = np.array((0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6))/60
ocean_color = np.float64([209,230,241])/255
land_color  = np.float64([0.6,0.6,0.6])

# Create New Colormap With a Limited Number of entries
nmap=plt.cm.get_cmap(name=plt.cm.Reds,lut=14) # only needed to set bins for pcolor/pcolormesh
nmap=plt.cm.get_cmap(name=plt.cm.OrRd,lut=14) # only needed to set bins for pcolor/pcolormesh

# Recenter over Pacific
lon_new=lon+180
i_west=np.where(lon<=0)[0]
i_east=np.where(lon>=0)[0]
sst_new = np.hstack((djf_sst_trend[:,i_east],djf_sst_trend[:,i_west]))
                    
# Create Figure                
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=-180))
m = ax.contourf(lon_new, lat, sst_new, clevs, transform=ccrs.PlateCarree(),cmap=nmap,extend="max")
plt.title('Warming Trend (K/yr), 1860-1869',fontsize=28)
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.LAND,color=land_color,zorder=1)
#ax.text(10,-60,'Total Warming (K), 2014-2100',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal")
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs[np.arange(0,14,2)])
cbar.ax.tick_params(labelsize=24)
#cbar.ax.set_xticks(np.arange(0.0,0.06+1))
#plot_rectangle(ax,-120,-90,20,40)
plot_rectangle(ax,360-90,360-80,-10,0)
plot_rectangle(ax,160,360-150,-5,5)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/djf_sst_map_trend_2015_2100.eps',format='eps')



### DJF SST Anomalies: period 2 (relative to period 1, which could include period 2)
# Setup some map stuff
clevs       = np.array((-0.5,-0.4,-0.3,-0.2,-0.1,0,.1,.2,.3,.4,.5))
clevs       = np.array((-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,.1,.2,.3,.4,.5,.6))
ocean_color = np.float64([209,230,241])/255
land_color  = np.float64([0.6,0.6,0.6])

# Create New Colormap With a Limited Number of entries
nmap=plt.cm.get_cmap(name=plt.cm.RdBu_r,lut=12) # only needed to set bins for pcolor/pcolormesh

# 20th Century
yrs_drght_20th = np.arange(1863,1866)
i_dec    = np.where(  (yrs>=np.min(yrs_drght_20th)-1) & (yrs<=np.max(yrs_drght_20th)-1) )[0]
i_janfeb = np.where(  (yrs>=np.min(yrs_drght_20th))   & (yrs<=np.max(yrs_drght_20th)) )[0]

# Calculate DJF SST Anomaly for these years
anomSST_drght20th = np.mean((anom_sst20th[i_dec,11,:,:]+anom_sst20th[i_janfeb,0,:,:]+anom_sst20th[i_janfeb,1,:,:])/3 \
                           , axis=0)

# Recenter over Pacific
lon_new=lon+180
i_west=np.where(lon<=0)[0]
i_east=np.where(lon>=0)[0]
sst_new = np.hstack((anomSST_drght20th[:,i_east],anomSST_drght20th[:,i_west]))
                    
# Create Figure                
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=-180))
m = ax.contourf(lon_new, lat, sst_new, clevs, transform=ccrs.PlateCarree(),cmap=nmap,extend="both")
ax.coastlines()
ax.set_global()
ax.add_feature(cartopy.feature.LAND,color=land_color,zorder=1)
ax.text(10,-60,'1863-1866',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal")
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plot_rectangle(ax,360-90,360-80,-10,0)
plot_rectangle(ax,160,360-150,-5,5)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/djf_sst_map_anom_1948_1957.eps',format='eps')

### Load continuous historical simulations (1860--1869)
# Calculate climatologies for a period (1860--1869)

# Load Data----------------------------------------------------------------------------------------------------------
sim_name  = 'CarbonSHPs'; 
ens_list = ['D01']

# Unique yrs
yrs = np.arange(1860,1869+1)

# Full year/month vectors
yrs_vect  = np.repeat(np.arange(1860,1869+1),12,axis=0)
mons_vect = np.ndarray.flatten(np.repeat([np.arange(1,12+1)],np.unique(yrs_vect).size,axis=0))

# Storage Arrays (climatology stuff)
clim_temp        = np.zeros((12,lat.size,lon.size))*np.nan
clim_prec        = np.zeros((12,lat.size,lon.size))*np.nan
clim_evap        = np.zeros((12,lat.size,lon.size))*np.nan
clim_irrig       = np.zeros((12,lat.size,lon.size))*np.nan
clim_irrig_gw    = np.zeros((12,lat.size,lon.size))*np.nan
clim_mean_smsurf = np.zeros((12,lat.size,lon.size))*np.nan
clim_sdev_smsurf = np.zeros((12,lat.size,lon.size))*np.nan
clim_mean_smroot = np.zeros((12,lat.size,lon.size))*np.nan
clim_sdev_smroot = np.zeros((12,lat.size,lon.size))*np.nan

# Loop Through Each Ensemble
for i_ens in enumerate(ens_list):
       
    # Load information from current ensemble member-------------------------------------------------------------------
    ens_name = ens_list[i_ens[0]]   
    print(ens_name)

    # Temperature (degrees C)-----------------------------------------------------------------------------------------
    
    # Historical
    ncfile     = netCDF4.Dataset(dir_name+'/ta_Amon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')     
    tsurf_all = ncfile.variables['ta'][:]
    ncfile.close
      
    # Precipitation (mm/day)------------------------------------------------------------------------------------------
    ncfile     = netCDF4.Dataset(dir_name+'/pr_Amon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')   
    prec_all = ncfile.variables['pr'][:]
    ncfile.close
    
    # Evaporation---------------------------------------------------------------------------------------------------
    ncfile     = netCDF4.Dataset(dir_name+'/evap_Lmon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')     
    evap_all        = ncfile.variables['evap'][:]
    ncfile.close    
    
    # Soil moisture---------------------------------------------------------------------------------------------------
    dir_name  = '/Users/puma/Soil_TerraE_Diag/processed';
    ncfile     = netCDF4.Dataset(dir_name+'/soilW_Lmon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')     
    soilW        = ncfile.variables['soilW'][:]
    smsurf_all = np.nansum(soilW[:,0:2,:,:],axis=1)  # top 2 layers, ~27cm
    smroot_all = np.nansum(soilW[:,0:4,:,:],axis=1)  # top 4 layers, ~1 meter
    ncfile.close

    # Total irrigation---------------------------------------------------------------------------------------------------
    dir_name  = '/Users/puma/Soil_TerraE_Diag/processed';
    ncfile     = netCDF4.Dataset(dir_name+'/irrig_w_Lmon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')     
    irrig_all        = ncfile.variables['irrig_w'][:]
    ncfile.close

    # Irrigation from fossil groundwater---------------------------------------------------------------------------------------------------
    dir_name  = '/Users/puma/Soil_TerraE_Diag/processed';
    ncfile     = netCDF4.Dataset(dir_name+'/irrig_gw_Lmon_ModelE_CarbonSHPs_'+ens_name+'_186001186912.nc')     
    irrigGW_all = ncfile.variables['irrig_gw'][:]
    ncfile.close

### Calculate averages for one period
ta_meanlatlon = np.mean(tsurf_all[:,:,:],axis=0)
prec_meanlatlon = np.mean(prec_all[:,:,:],axis=0) 
evap_meanlatlon = np.mean(evap_all[:,:,:],axis=0) 
irrig_meanlatlon = np.mean(irrig_all[:,:,:],axis=0) 
irrigGW_meanlatlon = np.mean(irrigGW_all[:,:,:],axis=0) 
smsurf_meanlatlon = np.mean(smsurf_all[:,:,:],axis=0)
smroot_meanlatlon = np.mean(smroot_all[:,:,:],axis=0) 

#---------------------------------------------------------------------------------------------------------------------
extent_lonlat = (-135, -65, 15, 65)
latlon_swusa = np.array([-122,-104,28,37])
latlon_texas = np.array([-104,-93,28,37])
latlon_seusa = np.array([-93,-75,28,37])

# Precipitation (mm/day), 1860-1869
data_map = copy.deepcopy(prec_meanlatlon)
clevs        = 5*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
nmap=plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(newcoast, linewidth=1, linestyle='-', zorder=2,edgecolor='k')
ax.text(-134,17,'Precip. (annual)',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
ax.text(-132,14,'1860-1869',transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
# Add a Rectangle and some annotation
plot_rectangle(ax,latlon_swusa[0],latlon_swusa[1],latlon_swusa[2],latlon_swusa[3])
plot_rectangle(ax,latlon_texas[0],latlon_texas[1],latlon_texas[2],latlon_texas[3])
plot_rectangle(ax,latlon_seusa[0],latlon_seusa[1],latlon_seusa[2],latlon_seusa[3])
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs)
cbar.ax.tick_params(labelsize=20)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/precip_1860-1869.eps',format='eps')

# Create Figure and set projection
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
#ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/precipGlob_1860-1869.eps',format='eps')


# Evaporation (mm/day), 1860-1869 (shoudl be all compoents for evapotranspiration)
data_map = copy.deepcopy(evap_meanlatlon)
clevs        = 5*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
nmap=plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(newcoast, linewidth=1, linestyle='-', zorder=2,edgecolor='k')
ax.text(-134,17,'Evap. (annual)',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
ax.text(-132,14,'1860-1869',transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
# Add a Rectangle and some annotation
plot_rectangle(ax,latlon_swusa[0],latlon_swusa[1],latlon_swusa[2],latlon_swusa[3])
plot_rectangle(ax,latlon_texas[0],latlon_texas[1],latlon_texas[2],latlon_texas[3])
plot_rectangle(ax,latlon_seusa[0],latlon_seusa[1],latlon_seusa[2],latlon_seusa[3])
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs)
cbar.ax.tick_params(labelsize=20)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/evap_1860-1869.eps',format='eps')

# Create Figure and set projection
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
#ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/evapGlob_1860-1869.eps',format='eps')


# Temperature (C), 1860-1869
data_map = copy.deepcopy(ta_meanlatlon)
clevs        = np.array((-40,-30,-20,-10,0,10,20,30,40))
#nmap=plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh
nmap=plt.cm.get_cmap(name='viridis',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh
#nmap=plt.cm.get_cmap(name=plt.cm.BrBG,lut=9) # only needed to set bins for pcolor/pcolormesh
#nmap=plt.cm.get_cmap(name='BrBG',lut=12)

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(newcoast, linewidth=1, linestyle='-', zorder=2,edgecolor='k')
ax.text(-134,17,'Temperature (annual)',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
ax.text(-132,14,'1860-1869',transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
# Add a Rectangle and some annotation
plot_rectangle(ax,latlon_swusa[0],latlon_swusa[1],latlon_swusa[2],latlon_swusa[3])
plot_rectangle(ax,latlon_texas[0],latlon_texas[1],latlon_texas[2],latlon_texas[3])
plot_rectangle(ax,latlon_seusa[0],latlon_seusa[1],latlon_seusa[2],latlon_seusa[3])
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs)
cbar.ax.tick_params(labelsize=20)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/temp_1860-1869.eps',format='eps')

# Create Figure and set projection
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
#ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/tempGlob_1860-1869.eps',format='eps')

# Evaporation (mm/day), 1860-1869
data_map = copy.deepcopy(prec_meanlatlon)
clevs        = 5*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
nmap=plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(newcoast, linewidth=1, linestyle='-', zorder=2,edgecolor='k')
ax.text(-134,17,'Precip. (annual)',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
ax.text(-132,14,'1860-1869',transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
# Add a Rectangle and some annotation
plot_rectangle(ax,latlon_swusa[0],latlon_swusa[1],latlon_swusa[2],latlon_swusa[3])
plot_rectangle(ax,latlon_texas[0],latlon_texas[1],latlon_texas[2],latlon_texas[3])
plot_rectangle(ax,latlon_seusa[0],latlon_seusa[1],latlon_seusa[2],latlon_seusa[3])
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs)
cbar.ax.tick_params(labelsize=20)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/precip_1860-1869.eps',format='eps')

# Create Figure and set projection
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
#ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/precipGlob_1860-1869.eps',format='eps')


# Soil Moisture, ~1 meter (mm), 1860-1869
data_map = copy.deepcopy(smroot_meanlatlon)
clevs        = 600*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
nmap=plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(newcoast, linewidth=1, linestyle='-', zorder=2,edgecolor='k')
ax.text(-134,17,'Soil Moist. ~1 m (annual)',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
ax.text(-132,14,'1860-1869',transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
# Add a Rectangle and some annotation
plot_rectangle(ax,latlon_swusa[0],latlon_swusa[1],latlon_swusa[2],latlon_swusa[3])
plot_rectangle(ax,latlon_texas[0],latlon_texas[1],latlon_texas[2],latlon_texas[3])
plot_rectangle(ax,latlon_seusa[0],latlon_seusa[1],latlon_seusa[2],latlon_seusa[3])
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs)
cbar.ax.tick_params(labelsize=20)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/SM1m_1860-1869.eps',format='eps')

# Create Figure and set projection
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
#ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/SM1m_Glob_1860-1869.eps',format='eps')

# Soil Moisture, surface # top 2 layers, ~270mm, 1860-1869
data_map = copy.deepcopy(smsurf_meanlatlon)
clevs        = 200*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
nmap=plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=-100, central_latitude=25, globe=None))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
#ax.gridlines(xlocs=np.arange(-180,190,10),ylocs=np.arange(-180,190,10))
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
ax.add_feature(newcoast, linewidth=1, linestyle='-', zorder=2,edgecolor='k')
ax.text(-134,17,'Soil Moist. ~270mm (annual)',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
ax.text(-132,14,'1860-1869',transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
        horizontalalignment='left', verticalalignment='center',)
# Add a Rectangle and some annotation
plot_rectangle(ax,latlon_swusa[0],latlon_swusa[1],latlon_swusa[2],latlon_swusa[3])
plot_rectangle(ax,latlon_texas[0],latlon_texas[1],latlon_texas[2],latlon_texas[3])
plot_rectangle(ax,latlon_seusa[0],latlon_seusa[1],latlon_seusa[2],latlon_seusa[3])
# Colorbar
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs)
cbar.ax.tick_params(labelsize=20)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/SM27cm_1860-1869.eps',format='eps')

# Create Figure and set projection
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
#m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='viridis',extend="max")
m = ax.pcolormesh(lon-1.25, lat-1, data_map, transform=ccrs.PlateCarree(), \
              cmap=nmap, vmin=np.min(clevs),vmax=np.max(clevs))
ax.coastlines()
ax.set_global()
#ax.set_extent(extent_lonlat, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
#ax.add_feature(states_provinces, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newcoast, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(newlake, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.LAND,color='w',zorder=0,edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=0,edgecolor='k')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='-', edgecolor='k')
ax.add_feature(cartopy.feature.OCEAN,color=ocean_color,zorder=1)
cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04)
cbar.ax.tick_params(labelsize=24)
plt.show()
fig.savefig('/Users/puma/Soil_TerraE_Diag/processed/figures/SM27cm_Glob_1860-1869.eps',format='eps')
