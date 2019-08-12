# -------------------------------------------------------------------------
#
#     #######    ########   ####     ####  
#    ##     ##      ##     ##   ##  ##   ##
#    ##             ##     ###      ###    
#    ##             ##       ###      ###  
#    ##     ####    ##         ###      ###
#    ##     ##      ##     ##   ##  ##   ##
#     #######    ########   #####    ##### 
#
# Name:        alter_prec_data
#
# Purpose:     Change the Precipitation value of netCDF4 files used in the CWat
#              model. The code can be changed to suite other modifications to
#              the netCDF4 files.  
#
# Author:      N. Woodruff, M.J.Puma

# Created:     16/07/2019
# -------------------------------------------------------------------------

import calendar
import netCDF4

# Percentage change to be implemented to the precipitation variable
percentPr = 0.9

# Open up netCDF file
# "r+" is read the data, alter it, and save it back to the original file.
with netCDF4.Dataset("0.9_pr_rhine.nc", "r+") as dset:
    
    # Set the data of precipitation to be a new variable
    # The notes after the OriginalData are the different sizes and meanings of the different parts of the 3 dimensional array
    OriginalData = dset['prec'][:,:,:]  #[ time , lat , lon ] [19358, 12, 14]
    
    numDaysInRegYear = 365                 # Initialize the data variable
    numDaysInLeapYear = 366

    totalDays = 0                       # Initialize variable
    
    # Loop through all the years the netCDF file contains (You have to know this number)
    for year in range(1961, 2014):
        if calendar.isleap(year) == True:
            for j in range(numDaysInLeapYear):
                
                # Summer time -> June-July-Aug 
                # Alter the original data, and save it to the file
                if 152 <= j <= 243: 
                    dset['prec'][totalDays,:,:] = percentPr * OriginalData[totalDays,:,:] 
                
                # increment total days, to keep track of which day you are altering in the file
                totalDays += 1
        
        else:
            for k in range(numDaysInRegYear):
                
                # Summer time -> June-July-August
                # Alter the original data, and save it to the file
                if 151 <= k <= 242:
                    dset['prec'][totalDays,:,:] = percentPr * OriginalData[totalDays,:,:] 
                
                # increment total days, to keep track of which day you are altering in the file
                totalDays += 1
