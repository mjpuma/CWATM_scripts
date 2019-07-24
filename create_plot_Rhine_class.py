#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on July 17 2019
Code to analyze and plot CWATM diagnostics processed by script 
"create_timeseries_CWATM.py"
@author:  M. Puma; N. Woodruff
'''

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

#------------------------------------------------------------------------------
#                           Start of the class
# The purpose of this class is to take all the CWatM Rhine Output monthly data
# and make difference plots. The code can also be altered to create regular plots
#------------------------------------------------------------------------------
class Rhine_Meteo_Var_Plots:
    def __init__(self):
        ## Directories
        self.dir_name  = 'C:/Users/cryst/Desktop/nasa/rhine30min/output';

        #----------------------------------------------------------------------
        # I'm pretty sure nothing else gets used here, but just in the off chance that it is used
        #----------------------------------------------------------------------
        # Month Vector
        self.mons     = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    
    #--------------------------------------------------------------------------
    # For obtaining lat and lon values
    #--------------------------------------------------------------------------
    def load_lat_lon(self):
            # Pull out time invariant variables (gridcell area, lat, lon)
            fname = self.dir_name + '/orig/monthly/Precipitation_monthavg.nc'
                
            # Open this netcdf file
            ncfile = netCDF4.Dataset(fname) # Load Dimension/Fixed Variables
            
            # Load variables
            self.lat      = ncfile.variables['lat'][:]         # latitude,  degrees
            self.lon      = ncfile.variables['lon'][:]         # longitude, degrees
            
            # Close the netcdf file
            ncfile.close 

    # For plotting a rectangle on the maps
    def plot_rectangle(self,ax, lonmin,lonmax,latmin,latmax):
        xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
        ys = [latmin,latmin,latmax,latmax,latmin]
        #ax.plot(xs, ys,latlon = True, color='k', linestyle='--', linewidth=3)
        ax.plot(xs,ys,color='k',linestyle='--',linewidth=3,transform=ccrs.PlateCarree())

    def read_orig_var(self):
        # Load Original Data----------------------------------------------------------------------------------------------------------
        sim_name  = 'Rhine River Precipitation Sensitivity Study'; 
        #for enumerating, uncomment next line
        #ens_list = ['D01']
        
        # Unique yrs
        yrs = np.arange(1990,2013+1)
        
        # Full year/month vectors
        yrs_vect  = np.repeat(np.arange(1990,2013+1),12,axis=0)
        mons_vect = np.ndarray.flatten(np.repeat([np.arange(1,12+1)],np.unique(yrs_vect).size,axis=0))
        
        # Storage Arrays (climatology stuff)
        clim_Tavg           = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_ETRef          = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_runoff         = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_totalET        = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_baseflow       = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_discharge      = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_Precipitation  = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_sum_gwRecharge = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        
        # Temperature (degrees C)-----------------------------------------------------------------------------------------
        
        # Historical
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/Tavg_monthavg.nc')     
        Tavg_all   = ncfile.variables['Tavg_monthavg'][:,:,:]
        ncfile.close
        
        # Evapotranspiration---------------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/ETRef_monthavg.nc')     
        ETRef_all  = ncfile.variables['ETRef_monthavg'][:,:,:]
        ncfile.close    

        # runoff (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/runoff_monthavg.nc')     
        runoff_all = ncfile.variables['runoff_monthavg'][:,:,:]
        ncfile.close
        
        # totalET (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/totalET_monthavg.nc')     
        totalET_all = ncfile.variables['totalET_monthavg'][:,:,:]
        ncfile.close

        # baseflow (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/baseflow_monthavg.nc')     
        baseflow_all = ncfile.variables['baseflow_monthavg'][:,:,:]
        ncfile.close

        # discharge (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/discharge_monthavg.nc')     
        discharge_all = ncfile.variables['discharge_monthavg'][:,:,:]
        ncfile.close
        
        # Precipitation (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/Precipitation_monthavg.nc')     
        Precipitation_all = ncfile.variables['Precipitation_monthavg'][:,:,:]
        ncfile.close
        
        # sum Ground Water Recharge (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/orig/monthly/sum_gwRecharge_monthavg.nc')     
        sum_gwRecharge_all = ncfile.variables['sum_gwRecharge_monthavg'][:,:,:]
        ncfile.close
        
        ### Calculate averages for one period
        Tavg_meanlatlon = np.mean(Tavg_all[:,:,:],axis=0)
        ETRef_meanlatlon = np.mean(ETRef_all[:,:,:],axis=0) 
        runoff_meanlatlon = np.mean(runoff_all[:,:,:],axis=0) 
        totalET_meanlatlon = np.mean(totalET_all[:,:,:],axis=0) 
        baseflow_meanlatlon = np.mean(baseflow_all[:,:,:],axis=0) 
        discharge_meanlatlon = np.mean(discharge_all[:,:,:],axis=0)
        Precipitation_meanlatlon = np.mean(Precipitation_all[:,:,:],axis=0)
        sum_gwRecharge_meanlatlon = np.mean(sum_gwRecharge_all[:,:,:], axis = 0)

        self.orig_meteo_var = np.array([Tavg_meanlatlon, ETRef_meanlatlon, runoff_meanlatlon, totalET_meanlatlon, baseflow_meanlatlon, discharge_meanlatlon, Precipitation_meanlatlon, sum_gwRecharge_meanlatlon])
        
    def read_new_var(self):
        # Load Data----------------------------------------------------------------------------------------------------------
        sim_name  = 'Rhine River Precipitation Sensitivity Study'; 
        #for enumerating, uncomment next line
        #ens_list = ['D01']
        
        # Unique yrs
        yrs = np.arange(1990,2013+1)
        
        # Full year/month vectors
        yrs_vect  = np.repeat(np.arange(1990,2013+1),12,axis=0)
        mons_vect = np.ndarray.flatten(np.repeat([np.arange(1,12+1)],np.unique(yrs_vect).size,axis=0))
        
        # Storage Arrays (climatology stuff)
        clim_Tavg           = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_ETRef          = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_runoff         = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_totalET        = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_baseflow       = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_discharge      = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_Precipitation  = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        clim_sum_gwRecharge = np.zeros((12,self.lat.size,self.lon.size))*np.nan
        
        # Temperature (degrees C)-----------------------------------------------------------------------------------------
        
        # Historical
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/Tavg_monthavg.nc')     
        Tavg_all   = ncfile.variables['Tavg_monthavg'][:,:,:]
        ncfile.close
        
        # Evapotranspiration---------------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/ETRef_monthavg.nc')     
        ETRef_all  = ncfile.variables['ETRef_monthavg'][:,:,:]
        ncfile.close    

        # runoff (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/runoff_monthavg.nc')     
        runoff_all = ncfile.variables['runoff_monthavg'][:,:,:]
        ncfile.close
        
        # totalET (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/totalET_monthavg.nc')     
        totalET_all = ncfile.variables['totalET_monthavg'][:,:,:]
        ncfile.close

        # baseflow (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/baseflow_monthavg.nc')     
        baseflow_all = ncfile.variables['baseflow_monthavg'][:,:,:]
        ncfile.close

        # discharge (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/discharge_monthavg.nc')     
        discharge_all = ncfile.variables['discharge_monthavg'][:,:,:]
        ncfile.close
        
        # Precipitation (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/Precipitation_monthavg.nc')     
        Precipitation_all = ncfile.variables['Precipitation_monthavg'][:,:,:]
        ncfile.close
        
        # sum Ground Water Recharge (units?)------------------------------------------------------------------------------------------
        ncfile     = netCDF4.Dataset(self.dir_name+'/x1.2_Prec/monthly/sum_gwRecharge_monthavg.nc')     
        sum_gwRecharge_all = ncfile.variables['sum_gwRecharge_monthavg'][:,:,:]
        ncfile.close
        
        ### Calculate averages for one period
        Tavg_meanlatlon = np.mean(Tavg_all[:,:,:],axis=0)
        ETRef_meanlatlon = np.mean(ETRef_all[:,:,:],axis=0) 
        runoff_meanlatlon = np.mean(runoff_all[:,:,:],axis=0) 
        totalET_meanlatlon = np.mean(totalET_all[:,:,:],axis=0) 
        baseflow_meanlatlon = np.mean(baseflow_all[:,:,:],axis=0) 
        discharge_meanlatlon = np.mean(discharge_all[:,:,:],axis=0)
        Precipitation_meanlatlon = np.mean(Precipitation_all[:,:,:],axis=0)
        sum_gwRecharge_meanlatlon = np.mean(sum_gwRecharge_all[:,:,:], axis = 0)

        self.new_meteo_var = np.array([Tavg_meanlatlon, ETRef_meanlatlon, runoff_meanlatlon, totalET_meanlatlon, baseflow_meanlatlon, discharge_meanlatlon, Precipitation_meanlatlon, sum_gwRecharge_meanlatlon])

    
    def plot(self, diffPlot = True, origPlot = True):
        #----------------------------------------------------------------------
        # Set Up for Plot
        #----------------------------------------------------------------------
        
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
        i_JJA    = np.where((self.mons>=6) & (self.mons<=8))[0]
        i_OND    = np.where((self.mons>=10) & (self.mons<=12))[0]
        i_JFM    = np.where((self.mons>=1) & (self.mons<=3))[0]
        i_AMJJAS = np.where((self.mons>=4) & (self.mons<=9))[0]

        # Setup some map stuff
        clevs       = np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
        ocean_color = np.float64([209,230,241])/255
        
        # Lat/Lon range for regional map
        extent_lonlat = (5, 12, 46, 52) # (min lon, max lon, min lat, max lat)
        
        # Create New Colormap With a Limited Number of entries
        nmap=plt.cm.get_cmap(name='viridis',lut=9) # only needed to set bins for pcolor/pcolormesh
        #nmap=plt.cm.get_cmap(name=plt.cm.BrBG,lut=9) # only needed to set bins for pcolor/pcolormesh
        #nmap=plt.cm.get_cmap(name='BrBG',lut=12)
        
        
        #----------------------------------------------------------------------
        # Boolean statement to determine plots plotted
        #----------------------------------------------------------------------
        
        if diffPlot == True:
            meteo_var = self.new_meteo_var - self.orig_meteo_var
        
        elif diffPlot == False and origPlot == True:    
            meteo_var = self.orig_meteo_var
        
        elif diffPlot == False and origPlot == False:
            meteo_var = self.new_meteo_var

        #-----------------------------------------------------------------------
        # Loop through all the arrays of the meteo variables, creating a plot for every one
        #-----------------------------------------------------------------------

        for i in range(len(meteo_var)):
            
            #---------------------------------------------------------------------------------------------------------------------
            extent_lonlat = (3, 12, 45, 53)
           
            # Precipitation (mm/day), 1860-1869
            data_map = copy.deepcopy(meteo_var[i])
            if i == 0:
                clevs = 20*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
            elif i == 1:
                clevs = 3e-3*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
            elif i == 2:
                clevs = 5e-3*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
            elif i == 3:
                clevs = 2e-3*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
            elif i == 4:
                clevs = 5e-4*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
            elif i == 5:
                clevs = 6e2*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
            elif i == 6:
                clevs = 5e-3*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
            elif i == 7:
                clevs = 7e-4*np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
                
            nmap=plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh
            
            fig = plt.figure(figsize=(12, 12))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=6, central_latitude=49, globe=None))
            #m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
            m = ax.pcolormesh(self.lon-1.25, self.lat-1, data_map, transform=ccrs.PlateCarree(), \
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
            
            if i == 0:
                ax.text(4,54,'Annual Temperature',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            elif i == 1:
                ax.text(4,54,'Annual Evapotranspiration',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            elif i == 2:
                ax.text(4,54,'Annual Runoff',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            elif i == 3:
                ax.text(4,54,'Annual total Evaporation',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            elif i == 4:
                ax.text(4,54,'Annual baseflow',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            elif i == 5:
                ax.text(4,54,'Annual discharge',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            elif i == 6:
                ax.text(4,54,'Annual Precipitation',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            elif i == 7:
                ax.text(4,54,'Annual total Ground Water Recharge',transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            
            ax.text(6.5,53.5,'1990-2013',transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
            # Add a Rectangle and some annotation
            #plot_rectangle(ax,latlon_swusa[0],latlon_swusa[1],latlon_swusa[2],latlon_swusa[3])
            #plot_rectangle(ax,latlon_texas[0],latlon_texas[1],latlon_texas[2],latlon_texas[3])
            #plot_rectangle(ax,latlon_seusa[0],latlon_seusa[1],latlon_seusa[2],latlon_seusa[3])
            
            # Colorbar
            cbar=plt.colorbar(m,orientation="horizontal",fraction=0.08,pad=0.04,ticks=clevs)
            cbar.ax.tick_params(labelsize=10)
            plt.show()
            
            if i == 0:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/Tavg.jpg',format='JPG')
            elif i == 1:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/ETRef.jpg',format='JPG')
            elif i == 2:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/runoff.jpg',format='JPG')
            elif i == 3:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/totalET.jpg',format='JPG')
            elif i == 4:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/baseflow.jpg',format='JPG')
            elif i == 5:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/discharge.jpg',format='JPG')
            elif i == 6:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/Precipitation.jpg',format='JPG')
            elif i == 7:
                fig.savefig('C:/Users/cryst/Desktop/nasa/rhine30min/output/Plots/sum_gwRecharge.jpg',format='JPG')
            
            
#------------------------------------------------------------------------------
# Initialize Class
#------------------------------------------------------------------------------

# Decide now wether you want to enumerate or not
#        # Loop Through Each Ensemble
#        for i_ens in enumerate(ens_list):
#            
#            # Load information from current ensemble member-------------------------------------------------------------------
#            ens_name = ens_list[i_ens[0]]   
#            print(ens_name)
        
Prec_Change = Rhine_Meteo_Var_Plots()
Prec_Change.load_lat_lon()
Prec_Change.read_orig_var()
Prec_Change.read_new_var()
Prec_Change.plot(diffPlot = False, origPlot = True)
