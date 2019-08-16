#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on July 17 2019
Code to analyze and plot CWATM diagnostics processed by script 
"create_timeseries_CWATM.py"
@author:  M. Puma; N. Woodruff
'''

# Import Modules and define functions - Everything you could possibly need is imported
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
    def __init__(self, folder_data_type, folder_data_change, units):

        # Month Vector
        self.mons     = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    
        # naming to help with handling the computation of all the different hydrological changes
        self.folder_data_change = folder_data_change
        self.folder_data_type   = folder_data_type
        
        # Determine which conversion factor will be used
        if folder_data_type == 'annual':
            self.unit_conv = 1000*365     #Units = mm/year
        elif folder_data_type == 'monthly':
            # Use this if you are averaging only certain months, and want daily averaged data in the region
            self.unit_conv = 1000
            
            # Use this if you are averaging over the months and want the units to be /month
            #self.unit_conv = 1000*12     # Units = mm/month
        else:
            self.unit_conv = 1000        # Units = mm/day
            
        # Some units do not have a '/time', and only need this conversion
        self.unit_conv_noTime = 1000     # Units = mm
        self.units = units
        
        ## Directories
        self.dir_name  = 'C:/Users/nwoodruf/Desktop/nasa/rhine30min/output'
        self.dir_name_2 = 'C:/Users/nwoodruf/Desktop/nasa/rhine30min/output/' + folder_data_change + '/' + folder_data_type;
        
    #--------------------------------------------------------------------------
    # For obtaining lat and lon values
    #--------------------------------------------------------------------------
    def load_lat_lon(self):
            # Pull out time invariant variables (gridcell area, lat, lon)
            fname = self.dir_name + '/orig/monthly/Precipitation_monthavg.nc'
            
            # Open this netcdf file
            with netCDF4.Dataset(fname) as ncfile:# Load Dimension/Fixed Variables
                # Load variables
                self.lat      = ncfile.variables['lat'][:]         # latitude,  degrees
                self.lon      = ncfile.variables['lon'][:]         # longitude, degrees
            
    def read_orig_var(self):
        # Load Original Data----------------------------------------------------------------------------------------------------------
        
        # Soil Moisture First Level (mm)-------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/sum_w1_annualavg.nc') as ncfile:
            sum_w1_all   = ncfile.variables['sum_w1_annualavg'][:,:,:]
            sum_w1_all   *= self.unit_conv_noTime
        
        # Soil Moisture Second Level (mm)-------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/sum_w2_annualavg.nc') as ncfile:     
            sum_w2_all   = ncfile.variables['sum_w2_annualavg'][:,:,:]
            sum_w2_all   *= self.unit_conv_noTime
            
        # Temperature (degrees C) Historical-----------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/Tavg_annualavg.nc') as ncfile:
            Tavg_all = ncfile.variables['Tavg_annualavg'][:,:,:]   # units (deg C)
        
        # Evapotranspiration---------------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/ETRef_annualavg.nc') as ncfile:
            ETRef_all  = ncfile.variables['ETRef_annualavg'][:,:,:]  # Units (m/day)
            ETRef_all *= self.unit_conv                                    # units converted to mm/year

        # runoff ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/runoff_annualavg.nc') as ncfile:
            runoff_all = ncfile.variables['runoff_annualavg'][:,:,:] # Units (m/day)
            runoff_all *= self.unit_conv                                   # Units converted to mm/year
        
        # totalET ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/totalET_annualavg.nc') as ncfile:     
            totalET_all = ncfile.variables['totalET_annualavg'][:,:,:]   # Units (m/day)
            totalET_all *= self.unit_conv                                      # units converted to mm/year   

        # baseflow ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/baseflow_annualavg.nc') as ncfile:    
            baseflow_all = ncfile.variables['baseflow_annualavg'][:,:,:] # Units (m/day)
            baseflow_all *= self.unit_conv                                     # Units converted to mm/year

        # discharge (m^3/year)------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/discharge_annualavg.nc') as ncfile:  
            discharge_all = ncfile.variables['discharge_annualavg'][:,:,:]   # Unites (m^3/month)
        
        # Precipitation ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/Precipitation_annualavg.nc') as ncfile:    
            Precipitation_all = ncfile.variables['Precipitation_annualavg'][:,:,:]   # Units (m/day)
            Precipitation_all *= self.unit_conv                                            # Units converted to mm/year
        
        # Ground Water Recharg (m/day)------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name+'/orig/annual/sum_gwRecharge_annualavg.nc') as ncfile:     
            sum_gwRecharge_all = ncfile.variables['sum_gwRecharge_annualavg'][:,:,:]   # Units (m/day)
            sum_gwRecharge_all *= self.unit_conv                                             # Units converted to mm/year
       
        #-----------------------------------------------------------------------
        # June July August Analysis code
        #------------------------------------------------------------------------
'''
        # Splice the data here - JJA
        totMon = 0
        extMon = 0
        for i in range(1981, 2011):
            numMon = 12 # same as self.mons but that is an array
            
            for j in range(numMon):
                
                # Summer time -> June-July-Aug 
                if 5 <= j <= 7: 
                    Tavg_JJA[extMon,:,:] = Tavg_all[totMon,:,:]
                    ETRef_JJA[extMon,:,:] = ETRef_all[totMon,:,:]
                    runoff_JJA[extMon,:,:] = runoff_all[totMon,:,:]
                    totalET_JJA[extMon,:,:] = totalET_all[totMon,:,:]
                    baseflow_JJA[extMon,:,:] = baseflow_all[totMon,:,:] 
                    discharge_JJA[extMon,:,:] = discharge_all[totMon,:,:]
                    Precipitation_JJA[extMon,:,:] = Precipitation_all[totMon,:,:]
                    sum_gwRecharge_JJA[extMon,:,:] = sum_gwRecharge_all[totMon,:,:]
                    sum_w1_JJA[extMon,:,:] = sum_w1_all[totMon,:,:]
                    sum_w2_JJA[extMon,:,:] = sum_w2_all[totMon,:,:]
                    
                    extMon += 1
                    
                totMon += 1
        
        ### Calculate averages for one period
        sum_w1_meanlatlon = np.mean(sum_w1_JJA[:,:,:],axis=0)
        sum_w2_meanlatlon = np.mean(sum_w2_JJA[:,:,:],axis=0)
        Tavg_meanlatlon = np.mean(Tavg_JJA[:,:,:],axis=0)
        ETRef_meanlatlon = np.mean(ETRef_JJA[:,:,:],axis=0) 
        runoff_meanlatlon = np.mean(runoff_JJA[:,:,:],axis=0) 
        totalET_meanlatlon = np.mean(totalET_JJA[:,:,:],axis=0) 
        baseflow_meanlatlon = np.mean(baseflow_JJA[:,:,:],axis=0) 
        discharge_meanlatlon = np.mean(discharge_JJA[:,:,:],axis=0)
        Precipitation_meanlatlon = np.mean(Precipitation_JJA[:,:,:],axis=0)
        sum_gwRecharge_meanlatlon = np.mean(sum_gwRecharge_JJA[:,:,:],axis=0)
        p_et = Precipitation_meanlatlon - totalET_meanlatlon   
'''        
        
        ### Calculate averages for one period
        sum_w1_meanlatlon = np.mean(sum_w1_all[:,:,:],axis=0)
        sum_w2_meanlatlon = np.mean(sum_w2_all[:,:,:],axis=0)
        Tavg_meanlatlon = np.mean(Tavg_all[:,:,:],axis=0)
        ETRef_meanlatlon = np.mean(ETRef_all[:,:,:],axis=0) 
        runoff_meanlatlon = np.mean(runoff_all[:,:,:],axis=0) 
        totalET_meanlatlon = np.mean(totalET_all[:,:,:],axis=0) 
        baseflow_meanlatlon = np.mean(baseflow_all[:,:,:],axis=0) 
        discharge_meanlatlon = np.mean(discharge_all[:,:,:],axis=0)
        Precipitation_meanlatlon = np.mean(Precipitation_all[:,:,:],axis=0)
        sum_gwRecharge_meanlatlon = np.mean(sum_gwRecharge_all[:,:,:],axis=0)
        p_e = Precipitation_meanlatlon - totalET_meanlatlon

        self.orig_meteo_var = np.array([Tavg_meanlatlon, ETRef_meanlatlon, 
                                        runoff_meanlatlon, totalET_meanlatlon, 
                                        baseflow_meanlatlon, discharge_meanlatlon, 
                                        Precipitation_meanlatlon, sum_gwRecharge_meanlatlon, 
                                        p_e, sum_w1_meanlatlon, sum_w2_meanlatlon])
        
    def read_new_var(self):
        # Load Data----------------------------------------------------------------------------------------------------------

        # Temperature (degrees C) Historical-----------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/Tavg_annualavg.nc') as ncfile:
            Tavg_all = ncfile.variables['Tavg_annualavg'][:,:,:]   # units (deg C)
        
        # Evapotranspiration---------------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/ETRef_annualavg.nc') as ncfile:
            ETRef_all  = ncfile.variables['ETRef_annualavg'][:,:,:]  # Units (m/day)
            ETRef_all *= self.unit_conv                                    # units converted to mm/year

        # runoff ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/runoff_annualavg.nc') as ncfile:
            runoff_all = ncfile.variables['runoff_annualavg'][:,:,:] # Units (m/day)
            runoff_all *= self.unit_conv                                   # Units converted to mm/year
        
        # totalET ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/totalET_annualavg.nc') as ncfile:     
            totalET_all = ncfile.variables['totalET_annualavg'][:,:,:]   # Units (m/day)
            totalET_all *= self.unit_conv                                      # units converted to mm/year   

        # baseflow ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/baseflow_annualavg.nc') as ncfile:    
            baseflow_all = ncfile.variables['baseflow_annualavg'][:,:,:] # Units (m/day)
            baseflow_all *= self.unit_conv                                     # Units converted to mm/year

        # discharge (m^3/s)------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/discharge_annualavg.nc') as ncfile:  
            discharge_all = ncfile.variables['discharge_annualavg'][:,:,:]   # Unites (m^3/month)
        
        # Precipitation ------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/Precipitation_annualavg.nc') as ncfile:    
            Precipitation_all = ncfile.variables['Precipitation_annualavg'][:,:,:]   # Units (m/day)
            Precipitation_all *= self.unit_conv                                            # Units converted to mm/year
        
        # Ground Water Recharg (m/day)------------------------------------------------------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/sum_gwRecharge_annualavg.nc') as ncfile:     
            sum_gwRecharge_all = ncfile.variables['sum_gwRecharge_annualavg'][:,:,:]   # Units (m/day)
            sum_gwRecharge_all *= self.unit_conv                                             # Units converted to mm/year
       
        # Soil Moisture First Level (mm)-------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/sum_w1_annualavg.nc') as ncfile:
            sum_w1_all   = ncfile.variables['sum_w1_annualavg'][:,:,:]
            sum_w1_all   *= self.unit_conv_noTime
        
        # Soil Moisture Second Level (mm)-------------------------------------------
        with netCDF4.Dataset(self.dir_name_2+'/sum_w2_annualavg.nc') as ncfile:     
            sum_w2_all   = ncfile.variables['sum_w2_annualavg'][:,:,:]
            sum_w2_all   *= self.unit_conv_noTime
               
        #-----------------------------------------------------------------------
        # June July August Analysis code
        #------------------------------------------------------------------------
'''
        # Splice the data here - JJA
        totMon = 0
        extMon = 0
        for i in range(1981, 2011):
            numMon = 12 # same as self.mons but that is an array
            
            for j in range(numMon):
                
                # Summer time -> June-July-Aug 
                if 5 <= j <= 7: 
                    Tavg_JJA[extMon,:,:] = Tavg_all[totMon,:,:]
                    ETRef_JJA[extMon,:,:] = ETRef_all[totMon,:,:]
                    runoff_JJA[extMon,:,:] = runoff_all[totMon,:,:]
                    totalET_JJA[extMon,:,:] = totalET_all[totMon,:,:]
                    baseflow_JJA[extMon,:,:] = baseflow_all[totMon,:,:] 
                    discharge_JJA[extMon,:,:] = discharge_all[totMon,:,:]
                    Precipitation_JJA[extMon,:,:] = Precipitation_all[totMon,:,:]
                    sum_gwRecharge_JJA[extMon,:,:] = sum_gwRecharge_all[totMon,:,:]
                    sum_w1_JJA[extMon,:,:] = sum_w1_all[totMon,:,:]
                    sum_w2_JJA[extMon,:,:] = sum_w2_all[totMon,:,:]
                    
                    extMon += 1
                    
                totMon += 1
        
        ### Calculate averages for one period
        sum_w1_meanlatlon = np.mean(sum_w1_JJA[:,:,:],axis=0)
        sum_w2_meanlatlon = np.mean(sum_w2_JJA[:,:,:],axis=0)
        Tavg_meanlatlon = np.mean(Tavg_JJA[:,:,:],axis=0)
        ETRef_meanlatlon = np.mean(ETRef_JJA[:,:,:],axis=0) 
        runoff_meanlatlon = np.mean(runoff_JJA[:,:,:],axis=0) 
        totalET_meanlatlon = np.mean(totalET_JJA[:,:,:],axis=0) 
        baseflow_meanlatlon = np.mean(baseflow_JJA[:,:,:],axis=0) 
        discharge_meanlatlon = np.mean(discharge_JJA[:,:,:],axis=0)
        Precipitation_meanlatlon = np.mean(Precipitation_JJA[:,:,:],axis=0)
        sum_gwRecharge_meanlatlon = np.mean(sum_gwRecharge_JJA[:,:,:],axis=0)
        p_et = Precipitation_meanlatlon - totalET_meanlatlon   
'''     
            
        ### Calculate averages for one period
        sum_w1_meanlatlon = np.mean(sum_w1_all[:,:,:],axis=0)
        sum_w2_meanlatlon = np.mean(sum_w2_all[:,:,:],axis=0)
        Tavg_meanlatlon = np.mean(Tavg_all[:,:,:],axis=0)
        ETRef_meanlatlon = np.mean(ETRef_all[:,:,:],axis=0) 
        runoff_meanlatlon = np.mean(runoff_all[:,:,:],axis=0) 
        totalET_meanlatlon = np.mean(totalET_all[:,:,:],axis=0) 
        baseflow_meanlatlon = np.mean(baseflow_all[:,:,:],axis=0) 
        discharge_meanlatlon = np.mean(discharge_all[:,:,:],axis=0)
        Precipitation_meanlatlon = np.mean(Precipitation_all[:,:,:],axis=0)
        sum_gwRecharge_meanlatlon = np.mean(sum_gwRecharge_all[:,:,:],axis=0)
        p_e = Precipitation_meanlatlon - totalET_meanlatlon

        self.new_meteo_var = np.array([Tavg_meanlatlon, ETRef_meanlatlon, 
                                       runoff_meanlatlon, totalET_meanlatlon, 
                                       baseflow_meanlatlon, discharge_meanlatlon, 
                                       Precipitation_meanlatlon, sum_gwRecharge_meanlatlon, 
                                       p_e, sum_w1_meanlatlon, sum_w2_meanlatlon])

    #-----------------------------------------------------------------------------------------------
    # Loop through all the arrays of the meteo variables, creating a plot for every one
    #---------------------------------------------------------------------------------------------------
    def scatter_plot(self):
        
        # 3 = totalET, 5 = discharge/streamflow, 6 = prec, 7 = GWRecharge
        # 
        meteo_var = np.array([self.new_meteo_var[3], self.new_meteo_var[5], 
                              self.new_meteo_var[6], self.new_meteo_var[7],
                              self.new_meteo_var[9], self.new_meteo_var[10]])

        meanE = np.mean(meteo_var[0])
        meanP = np.mean(meteo_var[2])
        meanP_E = meanP - meanE
        mean = meteo_var[1]
        meanDis = mean[0,2]
        meanGW = np.mean(meteo_var[3])
        meanSum1 = np.mean(meteo_var[4])
        meanSum2 = np.mean(meteo_var[5])
            
        return(np.array([meanP_E, meanDis, meanE, meanGW, meanSum1, meanSum2])) # P-E, streamflow, Evapotrans, Groundwater
    
    #----------------------------------------------------------------------------------------------------------
    # Function for creating geographical plots of mean yearly values.
    #-------------------------------------------------------------------------------------------------------------       
    def plot(self, diffPlot = True, origPlot = True):
        
        # Borders for mapping: Cultural borders
        states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='50m',
                facecolor='none')
        
        # Coastline
        newcoast = cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                                edgecolor='k', facecolor='none')
        #Lakes
        newlake = cfeature.NaturalEarthFeature('physical', 'lakes', '10m',
                                               edgecolor='k', facecolor='none')
        # Set plot styles
        # Formatting for titles
        fontdict_title = {'fontsize': 36}
        fig_size = np.array([10,10])
        
        # Formatting for figures
#        style_new = {'xtick.direction': 'in', \
#                     'ytick.direction': 'in', \
#                     'font.sans-serif': 'Arial'}
        
        # Regions for averaging (southwest & southeast)
#        latlon_swusa = np.array([-115,-100,28,40])
#        latlon_seusa = np.array([-95,-75,28,37])

        # Setup some map ocean color
        ocean_color = np.float64([209,230,241])/255
        clevsArray = np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
        
        # Create New Colormap With a Limited Number of entries
        nmap=plt.cm.get_cmap(name='viridis',lut=9) # only needed to set bins for pcolor/pcolormesh
        #nmap=plt.cm.get_cmap(name=plt.cm.BrBG,lut=9) # only needed to set bins for pcolor/pcolormesh
        #nmap=plt.cm.get_cmap(name='BrBG',lut=12)
        
        
        #----------------------------------------------------------------------
        # Boolean statement to determine plots plotted
        #----------------------------------------------------------------------
        
        if diffPlot == True:
            meteo_var = self.new_meteo_var - self.orig_meteo_var
            title = 'Difference '
        
        elif diffPlot == False and origPlot == True:    
            meteo_var = self.orig_meteo_var
            title = 'Original '
        
        elif diffPlot == False and origPlot == False:
            meteo_var = self.new_meteo_var
            title = 'New '

        #-----------------------------------------------------------------------
        # Loop through all the arrays of the meteo variables, creating a plot for every one
        #-----------------------------------------------------------------------
            
        # Lat/Lon range for regional map
        extent_lonlat = (-10, 25, 43, 58)
        
        for i in range(2,3):#len(meteo_var)): # the '-1' will make the sum gw recharge plot not print. 
            
            #---------------------------------------------------------------------------------------------------------------------
            # Hydrologic Variables 1990 - 2013
            data_map = copy.deepcopy(meteo_var[i])
            if i == 0: #Temp
                clevs = 12*clevsArray
            elif i == 1: #ETRef
                clevs = 8e2*clevsArray
            elif i == 2: #Runoff
                clevs = 2e-2*clevsArray
            elif i == 3: #ET
                clevs = 750*clevsArray
            elif i == 4: #Baseflow
                clevs = .04*clevsArray
            elif i == 5: #Discharge
                clevs = 5e2*clevsArray
            elif i == 6: #Prec
                clevs = 2e3*clevsArray
            elif i == 7: #Ground Water
                clevs = 1e2*clevsArray
            elif i == 8: #P-E
                clevs = 1e3*clevsArray
            elif i == 9: #50 mm > Moisture
                clevs = 25*clevsArray
            elif i == 10: # 2000 mm > Soil Moisture > 50 mm
                clevs = 350*clevsArray
    
            bluemap = plt.cm.get_cmap(name='Blues',lut=clevs.size-1) # only needed to set bins for pcolor/pcolormesh
            greenmap = plt.cm.get_cmap(name='Greens',lut=clevs.size-1)
            
            fig = plt.figure(figsize=(12, 12))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=6, central_latitude=49, globe=None))
            #m = ax.contourf(lon, lat, data_map, clevs, transform=ccrs.PlateCarree(),cmap='Blues',extend="max")
           
            m1 = ax.pcolormesh(self.lon-1.25, self.lat-1, data_map, transform=ccrs.PlateCarree(), \
                              cmap=bluemap, vmin=np.min(clevs),vmax=np.max(clevs))
            m2 = ax.pcolormesh(self.lon-1.25, self.lat-1, data_map, transform=ccrs.PlateCarree(), \
                              cmap= greenmap, vmin=np.min(clevs),vmax=np.max(clevs))
#            #9-class BuGn
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
                ax.text(-15,59, title + 'Annual Temperature ' + self.units[i], transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m1,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 1:
                ax.text(-15,59,title + 'Annual Ref. Evap. ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m2,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 2:
                ax.text(-15,59,title + self.folder_data_type + ' Runoff ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m1,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 3:
                ax.text(-15,59,title + self.folder_data_type + ' Baseline ET' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m2,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 4:
                ax.text(-15,59,title + self.folder_data_type + ' baseflow ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar)#
                cbar=plt.colorbar(m1,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 5:
                ax.text(-15,59, title + self.folder_data_type + ' Baseline Discharge ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m1,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 6:
                ax.text(-15,59, title + self.folder_data_type + ' Baseline Precipitation ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar)#
                cbar=plt.colorbar(m1,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 7:
                ax.text(-15,59.5, title + self.folder_data_type + ' Baseline Groundwater Recharge' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m1,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
                ax.text(-14.75,58.5,units,transform=ccrs.PlateCarree(),fontsize=28,fontweight="normal", \
                        horizontalalignment='left', verticalalignment='center',)

            elif i == 8:
                ax.text(-15,59,title + self.folder_data_type + '  Baseline P-E ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m2,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 9:
                ax.text(-15,59,title + self.folder_data_type + ' Soil Moisture 0-50 ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m2,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
            elif i == 10:
                ax.text(-15,59,title + self.folder_data_type + ' Soil Moisture 50-2000 ' + self.units[i],transform=ccrs.PlateCarree(),fontsize=32,fontweight="normal", \
                    horizontalalignment='left', verticalalignment='center',)
                # Colorbar
                cbar=plt.colorbar(m2,orientation="horizontal",fraction=0.08,pad=0.04)#,ticks=clevs)
                cbar.ax.tick_params(labelsize=15)
                
# Saving the plots as a JPG                
#            if i == 0:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'Tavg.jpg',format='JPG')
#            elif i == 1:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'ETRef.jpg',format='JPG')
#            elif i == 2:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'runoff.jpg',format='JPG')
#            elif i == 3:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'totalET.jpg',format='JPG')
#            elif i == 4:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'baseflow.jpg',format='JPG')
#            elif i == 5:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'discharge.jpg',format='JPG')
#            elif i == 6:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'Precipitation.jpg',format='JPG')
#            elif i == 7:
#                fig.savefig(self.dir_name + '/Plots/' + self.folder_data_change + '/' + title+'sum_gwRecharge.jpg',format='JPG')
            
#------------------------------------------------------------------------------
# Initialize Class
#------------------------------------------------------------------------------

# Decide now wether you want to enumerate or not
# If you choose not to enumerate through, then pick what variables you want to create the plots
P_ET = []
streamflow = []
Groundwater = []
ET = []
w1 = []
w2 = []

x = np.array([-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50])

folder_data_type = np.array(['annual', 'monthly', 'daily'])
folder_data_change = np.array(['x0.5_Prec', 'x0.6_Prec', 'x0.7_Prec', 'x0.8_Prec', 
                               'x0.9_Prec', 'x1.0_Prec', 'x1.1_Prec', 'x1.2_Prec',
                               'x1.3_Prec', 'x1.4_Prec', 'x1.5_Prec'])
#------------------------------------------------------------------------------
# Depending on the units you are plotting your graphs in, comment the lines of 
# code that you do not need.
#------------------------------------------------------------------------------
    
# Annual Units
units = np.array(['deg C', 'mm/year', 'mm/year', 'mm/year', 'm^3/s', 'mm/year',
                  'mm/year','mm/year','mm','mm'])
# Monthly Units
# units = np.array(['deg C', 'mm/month', 'mm/month', 'mm/month', 'm^3/s', 'mm/month',
#                  'mm/month','mm/month','mm','mm'])
# Daily Units    
# units = np.array(['deg C', 'mm/day', 'mm/day', 'mm/day', 'm^3/s', 'mm/day',
 #                 'mm/day','mm/day','mm','mm'])
    
    
#Prec_Change = Rhine_Meteo_Var_Plots(folder_data_type[0], folder_data_change[5], units)
#Prec_Change.load_lat_lon()
#Prec_Change.read_orig_var()
#Prec_Change.read_new_var()
#Prec_Change.plot(diffPlot = False, origPlot = True)


for i in range(len(folder_data_change)):
    Prec_Change = Rhine_Meteo_Var_Plots(folder_data_type[0], folder_data_change[i], units)
    Prec_Change.load_lat_lon()
    Prec_Change.read_orig_var()
    Prec_Change.read_new_var()
#    Prec_Change.plot(diffPlot = False, origPlot = True)
    data = Prec_Change.scatter_plot() # P-E, streamflow, Groundwater
    P_ET.append(data[0])
    streamflow.append(data[1])
    ET.append(data[2])
    Groundwater.append(data[3])
    w1.append(data[4])
    w2.append(data[5])
    
plt.figure(figsize=(8, 5), dpi=100)

plt.scatter(x, Groundwater)
plt.title('Annual Groundwater Recharge',fontsize=16,fontweight="normal")
plt.xlabel('Precipitation Difference (%)',fontsize=16,fontweight="normal")
plt.ylabel('mm/yr',fontsize=16,fontweight="normal")
plt.grid(True, which = 'both')
plt.axhline(y=0, color = 'k')
plt.axvline(x=0, color = 'k')
plt.show()

plt.figure(figsize=(8, 5), dpi=100)

plt.scatter(x, P_ET)
plt.title('Annual P - E',fontsize=16,fontweight="normal")
plt.xlabel('Precipitation Difference (%)',fontsize=16,fontweight="normal")
plt.ylabel('mm/yr',fontsize=16,fontweight="normal")
plt.grid(True, which = 'both')
plt.axhline(y=0, color = 'k')
plt.axvline(x=0, color = 'k')
plt.show()

plt.figure(figsize=(8, 5), dpi=100)

plt.scatter(x, streamflow)
plt.title('Annual Discharge at River Mouth',fontsize=16,fontweight="normal")
plt.xlabel('Precipitation Difference (%)',fontsize=16,fontweight="normal")
plt.ylabel('m^3/s',fontsize=16,fontweight="normal")
plt.grid(True, which = 'both')
plt.axhline(y=0, color = 'k')
plt.axvline(x=0, color = 'k')
plt.show()

plt.figure(figsize=(8, 5), dpi=100)

plt.scatter(x, ET)
plt.title('Annual Evapotranspiration',fontsize=16,fontweight="normal")
plt.xlabel('Precipitation Difference (%)',fontsize=16,fontweight="normal")
plt.ylabel('mm/yr',fontsize=16,fontweight="normal")
plt.grid(True, which = 'both')
#plt.axhline(y=0, color = 'k')
plt.axvline(x=0, color = 'k')
plt.show()

plt.figure(figsize=(8, 5), dpi=100)

plt.scatter(x, w1)
plt.title('Summer soil moisture (0-50 mm)',fontsize=16,fontweight="normal")
plt.xlabel('Precipitation Difference (%)',fontsize=16,fontweight="normal")
plt.ylabel('mm',fontsize=16,fontweight="normal")
plt.grid(True, which = 'both')
#plt.axhline(y=0, color = 'k')
plt.axvline(x=0, color = 'k')
plt.show()

plt.figure(figsize=(8, 5), dpi=100)

plt.scatter(x, w2)
plt.title('Summer soil moisture (50 - 2000 mm)',fontsize=16,fontweight="normal")
plt.xlabel('Precipitation Difference (%)',fontsize=16,fontweight="normal")
plt.ylabel('mm',fontsize=16,fontweight="normal")
plt.grid(True, which = 'both')
#plt.axhline(y=0, color = 'k')
plt.axvline(x=0, color = 'k')
plt.show()

#            if i == 0:
#                fig.savefig(self.dir_name + '/Plots/scatter/' + self.folder_data_change + '/' + title+'runoff_scatter.jpg',format='JPG')
#            elif i == 1:
#                fig.savefig(self.dir_name + '/Plots/scatter/' + self.folder_data_change + '/' + title+'totalET_scatter.jpg',format='JPG')
#            elif i == 2:
#                fig.savefig(self.dir_name + '/Plots/scatter/' + self.folder_data_change + '/' + title+'streamflow_scatter.jpg',format='JPG')
#            elif i == 3:
#                fig.savefig(self.dir_name + '/Plots/scatter/' + self.folder_data_change + '/' + title+'GW_scatter.jpg',format='JPG')
            
