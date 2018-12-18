#################################################
# esens_subsetting.py
#
# Contains methods to subset an ensemble
#
# Slice information
# ----------------
# smat/mmat (sensitivity vars and means)
#  0 - gph300
#  1 - gph500
#  2 - gph700
#  3 - gph850
#  4 - t300
#  5 - t500
#  6 - t700
#  7 - t850
#  8 - t925
#  9 - u300
#  10 - u500
#  11 - u700
#  12 - u850
#  13 - u925
#  14 - v300
#  15 - v500
#  16 - v700
#  17 - v850
#  18 - v925
#  19 - q850
#  20 - slp
#  21 - t2
#  22 - q2
#  23 - td2
#  24 - u10
#  25 - v10
#########################################################################



import numpy as np
from netCDF4 import Dataset
#import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#from cartopy import feature as cfeat
#import cmocean

def point(matrix):
    '''
    Builds a mask that is true for a matrix everywhere
    except the highest magnitude value point in the matrix.
    '''
    maxval = []
    mask = np.zeros_like(matrix, dtype=bool)
    matrix = np.abs(matrix)
    for i in range(len(matrix[:,0,0])):
        maxval.append(np.ma.max(matrix))
        point = np.where(matrix == maxval)
        mask[i] = (matrix != matrix[point])
    return mask, maxval

def percent(matrix, percent):
    '''
    Builds a mask that is true for the lowest
    n percent of values in a matrix.
    '''
    maxval = []
    mask = np.zeros_like(matrix, dtype=bool)
    matrix = np.abs(matrix)
    for i in range(len(matrix[:,0,0])):
        print(np.max(matrix[i]))
        threshval = np.percentile(matrix[i].compressed(), percent)
        print("Sens Threshold: ", percent)
        print("Threshold value: ", threshval)
        # Item to be ignored in calculation if it's greater than lower percentile
        #  and less than the upper percentile
        mask[i] = (matrix[i] < threshval)
        #print(np.where(~mask))
        maxval.append(np.max(matrix[i]))
        print("Max values: ", maxval)
    return mask, maxval

def ensSubset(wrfsensfile, analysis, memvalsfile, fullensnum,
              newensnum, sensvars=['500_hPa_GPH'], method=3, sens_thresh=70.):
    '''
    Uses the WRF sensitivity and interpolated analysis file as well as the
    sensitivity variable values for each member to determine the optimal
    ensemble subset given desired subset size and subsetting method.

    Inputs
    ------
    wrfsensfile --- String containing path to sensitivity NC file. Assumes
                    sensitivity data is stored in 'P_HYD' as dictated by
                    esensSPC.f.
    analysis ------ String containing path to interpolated analysis file.
    memvalsfile --- String containing path to output of sensvector.f that
                    stores sensitivity variable values for each member.
    fullensnum ---- Integer describing the number of members in the full
                    ensemble.
    newensnum ----- Integer describing subset size.
    method -------- Integer from 1-3 to choose subsetting technique.
                    (1) Automatically chooses the single highest sensitivity
                        grid point and subsets based on the lowest errors
                        at this grid point.
                    (2) Weights all (or some percentage of)
                        sensitivity variable fields by the
                        sensitivity field and subsets based on lowest
                        total error.
                    (3) Takes the top n% (optimal number tbd) of the
                        sensitivity field and only calculates and sums
                        errors at these grid points, subsetting on
                        lowest total error.
    Outputs
    -------
    Returns an integer list of subset ensemble members as well as a
    list of the names of sensitivity variables used.
    '''
    # Define varkeys for analysis and sensinds for WRF
    var = {'300_hPa_GPH': 0, '500_hPa_GPH': 1, '700_hPa_GPH': 2, '850_hPa_GPH': 3,
       '300_hPa_T': 4, '500_hPa_T': 5, '700_hPa_T': 6, '850_hPa_T': 7, '925_hPa_T': 8,
       '300_hPa_U-Wind': 9, '500_hPa_U-Wind': 10, '700_hPa_U-Wind': 11,
       '850_hPa_U-Wind': 12, '925_hPa_U-Wind': 13, '300_hPa_V-Wind': 14,
       '500_hPa_V-Wind': 15, '700_hPa_V-Wind': 16, '850_hPa_V-Wind': 17,
       '925_hPa_V-Wind': 18, '850_hPa_Q': 19, 'SLP': 20, '2m_Temp': 21,
       '2m_Q': 22, '2m_Dewpt': 23, '10m_U-Wind': 24, '10m_V-Wind': 25}
    sensinds = [var[sensvar] for sensvar in sensvars]
    print("Sensitivity indices:", sensinds)

    # Pull sensitivity field with dimensions (sensitivity variable index, ydim, xdim)
    wrfsens = Dataset(wrfsensfile)
    smat = wrfsens.variables['P_HYD']
    lons, lats = wrfsens.variables['XLONG'][0], wrfsens.variables['XLAT'][0]
    sensmat = smat[0,sensinds[:],:,:]
    # Mask of underground and missing data
    missing = (sensmat >= 9e9)
    #print("Max sens: ", np.ma.amax(sensmat_masked))
    # Pull analysis
    anl = Dataset(analysis)
    anlvar = np.ma.zeros((len(sensvars), len(lats[:,0]), len(lons[0,:])))
    for i in range(len(sensvars)):
        anlvar[i,:,:] = anl.variables[sensvars[i]][:,:]
    anl_missing = (anlvar >= 9e9)
    sensstrings = [sensvar.replace("_"," ") for sensvar in sensvars.copy()]
    sensmat_masked = np.ma.masked_array(sensmat, mask=(missing))

    try:
        print("Subset Method: ",str(method))
        if method == 1:
            # Mask by point
            mask, maxval = point(sensmat_masked)
        elif method == 2:
            # Mask by sensitivity threshold specific to proj method
            mask, maxval = percent(sensmat_masked, sens_thresh)
        elif method == 3:
            # Mask by percentage of field
            mask, maxval = percent(sensmat_masked, sens_thresh)
        else:
            raise NameError("Invalid method choice: {}".format(str(method)))
    except:
        dflt = 3
        print("Setting subsetting technique to {}".format(str(dflt)))
        mask, maxval = percent(sensmat_masked, sens_thresh)

#    for i in range(len(sensvars)):
#        print("Max sensitivity val for sens var {}: ".format(sensvars[i]),
#              np.ma.max(sens_masked[i]))
#        print("Min sensitivity val for sens var {}: ".format(sensvars[i]),
#              np.ma.min(sens_masked[i]))
#    #print(sensmat_masked[~sensmat_masked.mask])
#    #print(sens_masked[~sens_masked.mask])

    # New naming conventions for sens values netcdf
    varkeydict = {0: 'GPH_300', 1: 'GPH_500', 2: 'GPH_700', 3: 'GPH_850',
           4: 'T_300', 5: 'T_500', 6: 'T_700', 7: 'T_850', 8: 'T_925',
           9: 'U_300', 10: 'U_500', 11: 'U_700', 12: 'U_850', 13: 'U_925',
           14: 'V_300', 15: 'V_500', 16: 'V_700', 17: 'V_850', 18: 'V_925',
           20: 'SLP', 21: 'T2', 22: 'Q2', 23: 'TD2', 24: 'U10', 25: 'V10'}
    varkeys = [varkeydict[ind] for ind in sensinds]
    print("Sens keys for pulling from sensvector: ", varkeys)

    # Pull sensitivity variable values from each member
    memvals = Dataset(memvalsfile)
    memvar = np.zeros((len(varkeys), fullensnum, len(lats[:,0]), len(lons[0,:])))
    for k in range(len(varkeys)):
        memvar[k,:,:,:] = memvals.variables[varkeys[k]][:,:,:]

    # Start calculating differences between obs and members
    tmask = (missing | mask)
    #print(np.shape(tmask))
    error = np.ma.zeros((len(varkeys), fullensnum, len(lats[:,0]), len(lons[0,:])))
    # Apply percent sensfield and missing masks to analysis
    anlvar_masked = np.ma.masked_array(anlvar, mask=(tmask | anl_missing))
    #print("RAP value(s): ", anlvar_masked[tmask==False])
    sens_masked = np.ma.masked_array(sensmat_masked, mask=tmask)
    diff = np.zeros_like(error)
    for k in range(len(varkeys)):
        print("Min sensitivity val for sens var {}: ".format(sensvars[k]), np.ma.min(np.abs(sens_masked[k])))
        print("Max sensitivity val for sens var {}: ".format(sensvars[k]), np.ma.max(np.abs(sens_masked[k])))
        for i in range(fullensnum):
            # Apply percent sensfield and missing masks to member
            mem_masked = np.ma.masked_array(memvar[k,i], mask=(tmask[k]))
            if method == 2:
                #print("Max diff between mem and ob: ", np.ma.max(np.abs(mem_masked[:,:]-anlvar_masked[k,:,:])))
                diff[k,i,:,:] = mem_masked[:,:] - anlvar_masked[k,:,:]
                error[k,i,:,:] = np.ma.abs(sens_masked[k,:,:]*diff[k,i,:,:])
                #print("Max Error:", np.max(error[k,i,:,:]))
            else:
                error[k,i,:,:] = np.abs(mem_masked[:,:]-anlvar_masked[k,:,:])
            # Restructure for later mask
            error[k,i][tmask[k]] = np.NaN
    # Mask all error data that was masked from member fields
    error_masked = np.ma.masked_array(error, mask=(error == np.NaN))
    print("Min/Max absolute weighted error: ", np.nanmin(error_masked), np.nanmax(error_masked))
    print("Index of max abs weighted error: ", np.where(error_masked == np.nanmax(error_masked)))

    ####################################################
    # Test by plotting resulting sensitivity field
    ####################################################

#    # Get lat/lon data for plot extent
#    wrfrefpath = '/lustre/research/bancell/aucolema/HWT2016runs/2016051300/wrfoutREF'
#    wrfref = Dataset(wrfrefpath)
#    clon, clat = wrfref.CEN_LON, wrfref.CEN_LAT
#    tlat1, tlat2 = wrfref.TRUELAT1, wrfref.TRUELAT2
#
#    try:
#        # Plot
#        fig = plt.figure(figsize=(10, 10))
#
#        # Build projection/map
#        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=clon,
#                                                                               central_latitude=clat,
#                                                                               standard_parallels=(tlat1, tlat2)))
#        state_borders = cfeat.NaturalEarthFeature(category='cultural',
#                       name='admin_1_states_provinces_lakes', scale='50m', facecolor='None')
#        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
#        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
#        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
#        ax.set_extent([-120., -75., 20., 50.])
#        cflevs = np.linspace(-1*maxval[0], maxval[0], 101.)
#        #print(cflevs)
#        sensfield_threshold = ax.contourf(lons, lats, sensmat_masked[0,:,:],
#                                          vmin=cflevs[0], vmax=cflevs[-1], levels=cflevs, cmap='bwr',
#                                          transform=ccrs.PlateCarree())
#        fig.colorbar(sensfield_threshold, fraction=0.046, pad=0.04, ax=ax, orientation='horizontal',
#                         label='Sensitivity Magnitude')
#        plt.title('Sensitivity Field of ' + varkeys[0].replace("_"," "))
#        plt.savefig('sensfield.png')
#        plt.close()
#    except:
#        plt.close()
#        print('Plotting sensitivity field as a check failed.')
#
#    try:
#        # Plot
#        fig = plt.figure(figsize=(10, 10))
#
#        # Build projection/map
#        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=clon,
#                                                                               central_latitude=clat,
#                                                                               standard_parallels=(tlat1, tlat2)))
#        state_borders = cfeat.NaturalEarthFeature(category='cultural',
#                       name='admin_1_states_provinces_lakes', scale='50m', facecolor='None')
#        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
#        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
#        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
#        ax.set_extent([-120., -75., 20., 50.])
#        cflevs = np.linspace(-1*maxval[0], maxval[0], 101.)
#        #print(cflevs)
#        sensfield_threshold = ax.contourf(lons, lats, sens_masked[0,:,:],
#                                          vmin=cflevs[0], vmax=cflevs[-1], levels=cflevs, cmap='bwr',
#                                          transform=ccrs.PlateCarree())
#        fig.colorbar(sensfield_threshold, fraction=0.046, pad=0.04, ax=ax, orientation='horizontal',
#                         label='Sensitivity Magnitude')
#        plt.title('Masked Sensitivity Field of ' + varkeys[0].replace("_"," "))
#        plt.savefig('sensfield_threshold.png')
#        plt.close()
#    except:
#        plt.close()
#        print('Plotting masked sensitivity field as a check failed.')
## Onto Error fields
#
#    maxweighted = np.nanmax(np.absolute(error_masked[0,4,:,:]))
##
##    try:
##        # Now plot error field
##        fig = plt.figure(figsize=(10,10))
##        # Build projection/map
##        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
##        state_borders = cfeat.NaturalEarthFeature(category='cultural',
##                       name='admin_1_states_provinces_lakes', scale='50m', facecolor='None')
##        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
##        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
##        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
##        ax.set_extent([-120., -75., 20., 50.])
##        errormax = np.ma.amax(np.absolute(diff[0,4,:,:]))
##        print("Absolute Max for Member 5: ", errormax)
##        cflevs = np.linspace(0.0, maxweighted, 100.)
##        #print(cflevs)
##        #print(np.shape(error_masked))
##        errorfield = ax.contourf(lons, lats, np.absolute(diff[0,4,:,:]), cflevs,
##                                 cmap=cmocean.cm.amp, transform=ccrs.PlateCarree())
##        fig.colorbar(errorfield, fraction=0.046, pad=0.04, ax=ax, orientation='horizontal',
##                         label='Error')
##        plt.title('Error Field for Member 5 of ' + varkeys[0].replace("_"," "))
##        plt.savefig('unweighted_errorfield.png')
##        plt.close()
##    except:
##        plt.close()
##        raise
##        print('Plotting masked error field as a check failed.')
#
#    try:
#        # Now plot error field
#        fig = plt.figure(figsize=(10,10))
#        # Build projection/map
#        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
#        state_borders = cfeat.NaturalEarthFeature(category='cultural',
#                       name='admin_1_states_provinces_lakes', scale='50m', facecolor='None')
#        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
#        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
#        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
#        ax.set_extent([-120., -75., 20., 50.])
#        print("Absolute Max Weighted for Member 5: ", maxweighted)
#        cflevs = np.linspace(0.0, maxweighted, 100.)
#        #print(cflevs)
#        #print(np.shape(error_masked))
#        errorfield = ax.contourf(lons, lats, error_masked[0,4,:,:], cflevs, cmap=cmocean.cm.amp,
#                                     transform=ccrs.PlateCarree())
#        #errorfield.cmap.set_under(color=u'navy')
#        fig.colorbar(errorfield, fraction=0.046, pad=0.04, ax=ax, orientation='horizontal',
#                         label='Error')
#        if method == 2:
#            plt.title('Masked and Weighted Error Field for Member 5 of ' + varkeys[0].replace("_"," "))
#            plt.savefig('proj_errorfield.png')
#        else:
#            plt.title('Masked Error Field for Member 5 of ' + varkeys[0].replace("_"," "))
#            plt.savefig('rms_errorfield.png')
#        plt.close()
#    except:
#        plt.close()
#        raise
#        print('Plotting masked error field as a check failed.')

    # Sum total error and choose members with least error for subset
    summed_error = np.zeros((fullensnum))
    for i in range(fullensnum):
        summed_error[i] = np.nansum(error_masked[:,i,:,:])/np.ma.size(error_masked[:,i,:,:])
    print("Npts for member {}".format(1), np.ma.size(error_masked[:,0,:,:]))    
    print("Npts for member {}".format(i+1), np.ma.size(error_masked[:,i,:,:]))
    sorted_inds = summed_error.argsort()
    subset_mems = sorted_inds[:newensnum]+1 # Add one to correct zero-based
    print('Total errors: ', summed_error[sorted_inds])
    print('Subset mems: ', subset_mems)

    return subset_mems, sensstrings
