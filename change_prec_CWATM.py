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
# Name:        change_prec_CWATM
# Purpose:     Change the Precipitation value to study the effects of precipitation
#              on the surrounding area in the Rhine River. The code can be 
#              changed to suite other modifications to the netCDF4 files.  
#
# Author:      M. Puma; N. Woodruff
#
# Created:     16/07/2019
# -------------------------------------------------------------------------
import netCDF4 as nc

with nc.Dataset("pr_rhine_original.nc") as orig, nc.Dataset("pr_rhine.nc", "w") as new:
    
    # copy global attributes all at once via dictionary
    new.setncatts(orig.__dict__)
    
    # copy dimensions
    for name, dimension in orig.dimensions.items():
        new.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
        
    i = 0
                    
    # copy all file data except for the excluded
    for name, variable in orig.variables.items():
        i += 1
        new.createVariable(name, variable.datatype, variable.dimensions)
                
        if i == 4: #prec is the 4th variable printed
            new[name][:] = 1.2*orig[name][:]
        else:
            new[name][:] = orig[name][:]
                
#        if name == prec:
#            # copy variable attributes all at once via dictionary
#            new[name].setncatts(orig[name].__dict__)
