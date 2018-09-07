#!/bin/python3

import netCDF4 as nc
import numpy as np

dat1 = nc.Dataset('fss.nc')
dat2 = nc.Dataset('lastfew_fss.nc')

# Copy over 1st dataset
with dat1 as src, nc.Dataset('fssfinal.nc', "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)
print('Done w first dataset.')

# Copy over 2nd dataset
with dat2 as src, nc.Dataset('fssfinal.nc', "a") as dst:
     #print(dst)
     n = len(dst.variables['Run_Init'][:])
     for name, variable in src.variables.items():
         #print(name, variable)
         var = dst.variables[name]
         if len(np.shape(var)) > 1:
             dst[name][n:,:] = src[name][:]
         else:
             dst[name][n:] = src[name][:]
print('Done')
