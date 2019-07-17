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
# Name:        change prec CWATM
# Purpose:     Change the Precipitation value to study the effects of precipitation
#              on the surrounding area in the Rhine River. The code can be 
#              changed to suite other modifications to the netCDF4 files.  
#
# Author:      Nonnie Woodruff, intern at GISS 
# Mentor:      Dr. Michael J.Puma & James Miller

# Created:     16/07/2019
# -------------------------------------------------------------------------

from netCDF4 import Dataset
import numpy as np

path = "pr_rhine_original.nc"                    # The netCDF4 file in the same folder as this python script

with Dataset(path) as data:
    
    lon = data["lon"][:]
    lat = data["lat"][:]
    time = data["time"][:]
    prec = data["prec"][:,:,:]
    print(data)
    
path2 = "pr_rhine.nc"
data2 = Dataset(path2, 'w', format='NETCDF4_CLASSIC')

lonDim = data2.createDimension('lon',14)
latDim = data2.createDimension('lat',12)
timeDim = data2.createDimension('time',19358)

lonVar = data2.createVariable('lon', np.float32, dimensions = ("lon"))
latVar = data2.createVariable('lat', np.float32, dimensions = ("lat"))
timeVar = data2.createVariable('time', np.float32, dimensions = ("time"))
precVar = data2.createVariable('prec', np.float32, dimensions = ("time","lat","lon"))

lonVar[:] = lon[:]
latVar[:] = lat[:]
timeVar[:] = time[:]

# The meat of this program. Changing the precipitation before saving it to the .nc file
precVar[:,:,:] = 1.2*prec[:,:,:]

data2.close()