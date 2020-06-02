#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:49:55 2018

Compilation of verification statistic
calculations.

@author: aucolema
"""

import wrf
import numpy as np
from netCDF4 import Dataset, date2num
from datetime import timedelta
import xarray as xr
from datetime import datetime, timedelta
from subprocess import call
from scipy import ndimage
import pyproj
import scipy as sp
import os
import csv
import sys
from wrf_ens_tools.post import post_process as pp
# from profilehooks import profile

package_dir =  os.path.dirname(os.path.abspath(__file__))

def bilinear_interp(grid1x, grid1y, grid2x, grid2y, z):
    """
    A method which interpolates a function
    z(grid1x, grid1y) of a grid (grid1x, grid1y) to another
    grid (grid2x, grid2y). Returns an array from the approximated
    function of the second grid (approximation of z(grid2x, grid2y)).
    """
    # Pair flattened x and y values as coordinates
    coords_from = list(zip(grid1y.flatten(), grid1x.flatten()))
    Z = z.flatten()
    # Set up interpolation function with original grid and interp variable
    interp = sp.interpolate.LinearNDInterpolator(coords_from, Z, fill_value=9e9)
    # Interpolate to new grid
    interpolated_z = interp(grid2y, grid2x)

    return interpolated_z


def nearest_neighbor_ttu(runinitdate, sixhr, rtime, nbrhd=0.,
                         wrfrefpath='/lustre/research/bancell/aucolema/HWT2016runs'
                                    '/2016050800/wrfoutREFd2'):
    """
    Interpolates storm reports valid over a
    1-hr or 6-hr time frame to the
    native TTU WRF grid. Returns the WRF grid
    in the form of binary hits and misses based
    on SPC storm report locations.
    """
    # Get initialization date
    rdate = runinitdate + timedelta(hours=rtime)
    # Since reports are from 12Z - 1159Z, make sure
    #  we are grabbing the correct date.
    if (rtime > 23) & (runinitdate.hour == 12):
        runinitdatef = (runinitdate + timedelta(days=1)).strftime('%y%m%d')
    elif (rtime > 35) & (runinitdate.hour == 0):
        runinitdatef = (runinitdate + timedelta(days=1)).strftime('%y%m%d')
    else:
        runinitdatef = runinitdate.strftime('%y%m%d')
    print('Pullng SPC reports from ', runinitdatef)
    print('Response time ', rdate)

    # Using Robert Hepper's code for nearest neighbor
    # w/out calculating practically perfect
    # Get reports CSV file from web
    rptfile = runinitdatef+'_rpts_filtered.csv'
    add = 'www.spc.noaa.gov/climo/reports/'+rptfile
    call(['wget',add])

    # Make lists of report lats and lons
    try:
        with open(rptfile) as csvf:
            r = csv.reader(csvf)
            mylist = list(r)
    except IOError:
        print('Report CSV file could not be opened.')
        sys.exit()

    length = len(mylist)-3
    time = [0]*length
    lats = [0]*length
    lons = [0]*length
    ct = 0
    for f in mylist:
        if 'Time' not in f and 'Comments' not in f:
            time[ct] = int(str(f[0])[:2])
            lats[ct] = float(f[5])
            lons[ct] = float(f[6])
            ct = ct+1

    # Get WRF lats/lons as pperf grid
    dat = Dataset(wrfrefpath)
    lon = dat.variables['XLONG'][0]
    lat = dat.variables['XLAT'][0]
    # Convert DX to kilometers
    dx = dat.DX / 1000.
    dat.close()

    # If there aren't any reports, zero across grid
    if length == 0:
        grid = np.zeros_like(lon)

    # Otherwise, let's grid the reports
    else:
        # If six hour, mask reports by valid times in window
        if sixhr:
            hour = rdate.hour
            hours = [(hour - i)%24 for i in range(1,7)]
            mask = [(hr in hours) for hr in time]
            print('Reports valid {} to {}'.format(hours[-1],hours[0]))
        else:
            hour = rdate.hour
            mask = [(hr == (hour-1)%24) for hr in time]
            print('Reports valid {} to {}'.format((hour-1)%24,hour))
        try:
            # Set up empty grid onto correct projection
            grid = np.zeros_like(lon)
            NDFD = pyproj.Proj("+proj=lcc +lat_1=25 +lat_2=25 +lon_0=-95 +x_0=0 +y_0=0 +"
                               "ellps=WGS84 +datum=WGS84 +units=m +no_defs")

            # Convert lat/lon grid into projection space
            X, Y = NDFD(lon, lat)
            # Convert lat/lon reports into projection space
            x, y = NDFD(lons, lats)

            # Create KD-Tree for effecient lookup
            gpoints = np.array(list(zip(X.ravel(), Y.ravel())))
            gtree = sp.spatial.cKDTree(gpoints)

            # Run the KD-Tree to get distances to nearest gridbox and then
            # index of nearest grid point
            dists, inds = gtree.query(np.array(list(zip(x, y))),
                                      distance_upper_bound=1000000000.)

            # Convert index of 1D array into index of 2D lat/lon array
            xind, yind = np.unravel_index(inds[mask], X.shape)

            # Loop through all points and increment that grid cell by 1
            gridinds = np.indices(grid.shape)
            for xi, yi in zip(xind, yind):
                dists = np.sqrt(((gridinds[0,:,:]-xi)*dx)**2 + ((gridinds[1,:,:]-yi)*dx)**2)
                inds = np.where(dists <= nbrhd)
                grid[inds] = 1
        except:
            grid = np.zeros_like(lon)

        return grid

def calc_prac_perf(runinitdate, sixhr, rtime, nbrhd=0., sigma=2):
    """
    Implementation of SPC practically perfect
    calculations adapted from SPC code

    Inputs
    ------
    runinitdate -- datetime obj for
                    model initialization being
                    used.
    sixhr -------- boolean specifying whether
                    to calculate practically perfect
                    probs over six hr time window or
                    use one hr time window.
    rtime -------- time (in num fcst hrs
                    from runinit) to obtain six hr
                    or one hr practically perfect.
    nbrhd -------- optional float depicting searching
                    distance for LSRs on observation
                    grid
    sigma -------- optional float specifying sigma
                    to use for Gaussian filter.

    Outputs
    -------
    returns tuple containing pract perf probs
    and lons/lats (respectively) of pperf grid
    """
    #Get initialization date
    runinitdatef = runinitdate.strftime('%y%m%d')
    rdate = runinitdate + timedelta(hours=rtime)
    # Since reports are from 12Z - 1159Z, make sure
    #  we are grabbing the correct date.
    if (rtime > 23) & (runinitdate.hour == 12):
        runinitdatef = (runinitdate + timedelta(days=1)).strftime('%y%m%d')
    elif (rtime > 35) & (runinitdate.hour == 0):
        runinitdatef = (runinitdate + timedelta(days=1)).strftime('%y%m%d')
    print('Pullng SPC reports from ', runinitdatef)
    print('Response time ', rdate)

    #Get report CSV file from web
    rptfile = runinitdatef+'_rpts_filtered.csv'
    add = 'www.spc.noaa.gov/climo/reports/'+rptfile
    call(['wget',add])

    #Make lists of report lats and lons
    try:
        with open(rptfile) as csvf:
            r = csv.reader(csvf)
            mylist = list(r)
    except IOError:
        print('Report CSV file could not be opened.')
        sys.exit()

    length = len(mylist)-3
    time = [0]*length
    lats = [0]*length
    lons = [0]*length
    ct = 0
    for f in mylist:
        if 'Time' not in f and 'Comments' not in f:
            time[ct] = int(str(f[0])[:2])
            lats[ct] = float(f[5])
            lons[ct] = float(f[6])
            ct = ct+1

    # Get WRF lats/lons as pperf grid
    wrffile = '/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREFd2'
    dat = Dataset(wrffile)
    wrflon = dat.variables['XLONG'][0]
    wrflat = dat.variables['XLAT'][0]
    # Convert DX to kilometers
    dx = dat.DX / 1000.
    dat.close()

    #If there aren't any reports, practically perfect is zero across grid
    if length == 0:
        pperf = np.zeros_like(wrflon)

    #Otherwise, let's grid the reports
    else:
        # If six hour, mask reports by valid times in window
        if sixhr:
            hour = rdate.hour
            hours = [(hour - i)%24 for i in range(1,7)]
            mask = [(hr in hours) for hr in time]
            print('Reports valid {} to {}'.format(hours[-1], hours[0]))
        else:
            hour = rdate.hour
            mask = [(hr == (hour-1)%24) for hr in time]
            print('Reports valid {} to {}'.format((hour-1)%24,hour))
        try:
            #Set up empty grid onto correct projection
            grid = np.zeros_like(wrflon)
            NDFD = pyproj.Proj("+proj=lcc +lat_1=25 +lat_2=25 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

            # Convert lat/lon grid into projection space
            X, Y = NDFD(wrflon, wrflat)
            WRFX, WRFY = NDFD(wrflon, wrflat)

            # Convert lat/lon reports into projection space
            x, y = NDFD(lons, lats)

            #Create KD-Tree for effecient lookup
            gpoints = np.array(list(zip(X.ravel(), Y.ravel())))
            gtree = sp.spatial.cKDTree(gpoints)

            # Run the KD-Tree to get distances to nearest gridbox and then index of nearest grid point
            dists, inds = gtree.query(np.array(list(zip(x, y))), distance_upper_bound=1000000000.)

            # Convert index of 1D array into index of 2D lat/lon array
            xind, yind = np.unravel_index(inds[mask], X.shape)

            # Loop through all points and increment that grid cell by 1
            gridinds = np.indices(grid.shape)
            for xi, yi in zip(xind, yind):
                dists = np.sqrt(((gridinds[0,:,:]-xi)*dx)**2 + ((gridinds[1,:,:]-yi)*dx)**2)
                inds = np.where(dists <= nbrhd)
                grid[inds] = 1

            # Gaussian smoother over our grid to create practically perfect probs
            pperf = ndimage.gaussian_filter(grid, sigma=sigma, order=0)

        except:
            pperf = np.zeros_like(wrflon)

    #Remove report CSV file
    os.remove(rptfile)
    print("Practically Perfect min/max: ", np.min(pperf), ' , ', np.max(pperf))

    return pperf, wrflon, wrflat

def calc_prac_perf_spc_grid(runinitdate, sixhr, rtime, sigma=2,
                            wrfrefpath='/lustre/scratch/aucolema/2016052600/wrfoutREFd2'):
    """
    Implementation of practically perfect probability
    calculations adapted from Robert Hepper's code.
    Practically perfect probs are calculated on a grid
    with 80-km grid spacing and then interpolated to the
    native WRF grid.

    Inputs
    ------
    runinitdate --- datetime obj for
                    model initialization being
                    used.
    sixhr --------- boolean specifying whether
                    to calculate practically perfect
                    probs over six hr time window or
                    use one hr time window.
    rtime --------- time (in num fcst hrs
                    from runinit) to obtain six hr
                    or one hr practically perfect.
    sigma --------- optional float specifying sigma
                    to use for Gaussian filter.

    Outputs
    -------
    returns tuple containing pract perf probs
    and lons/lats (respectively) of pperf grid
    """
    #Get initialization date
    runinitdatef = runinitdate.strftime('%y%m%d')
    rdate = runinitdate + timedelta(hours=rtime)
    # Since reports are from 12Z - 1159Z, make sure
    #  we are grabbing the correct date.
    if (rtime > 23) & (runinitdate.hour == 12):
        runinitdatef = (runinitdate + timedelta(days=1)).strftime('%y%m%d')
    elif (rtime > 35) & (runinitdate.hour == 0):
        runinitdatef = (runinitdate + timedelta(days=1)).strftime('%y%m%d')
    print('Pullng SPC reports from ', runinitdatef)
    print('Response time ', rdate)

    # Get yesterday's reports CSV file from web
    rptfile = runinitdatef+'_rpts_filtered.csv'
    add = 'www.spc.noaa.gov/climo/reports/'+rptfile
    call(['wget',add])

    # Make lists of report lats and lons
    try:
        with open(rptfile) as csvf:
            r = csv.reader(csvf)
            mylist = list(r)
    except IOError:
        print('Report CSV file could not be opened.')
        sys.exit()

    length = len(mylist)-3
    time = [0]*length
    lats = [0]*length
    lons = [0]*length
    ct = 0
    for f in mylist:
        if 'Time' not in f and 'Comments' not in f:
            time[ct] = int(str(f[0])[:2])
            lats[ct] = float(f[5])
            lons[ct] = float(f[6])
            ct = ct+1

    # Get lats and lons for practically perfect grid
    ppfile = '{}/pperf_grid_template.npz'.format(package_dir)
    f = np.load(ppfile)
    lon = f["lon"]
    lat = f["lat"]
    f.close()

    # Get WRF lats/lons as pperf grid
    dat = Dataset(wrfrefpath)
    wrflon = dat.variables['XLONG'][0]
    wrflat = dat.variables['XLAT'][0]
    dat.close()

    # If there aren't any reports, practically perfect is zero across grid
    if length == 0:
        pperf = np.zeros_like(wrflon)
    # Otherwise, let's grid the reports
    else:
        # If six hour, mask reports by valid times in window
        if sixhr:
            hour = rdate.hour
            hours = [(hour - i) % 24 for i in range(1, 7)]
            mask = [(hr in hours) for hr in time]
            print("Reports valid {} to {}".format(hours[-1], hours[0]))
        else:
            hour = rdate.hour
            mask = [(hr == (hour-1)%24) for hr in time]
            print('Reports valid {} to {}'.format((hour-1) % 24, hour))
            # print('Report hours from reports used:', np.array(time)[mask])
            # print(mask)
            # print(time)
        try:
            # Set up empty grid onto correct projection
            grid = np.zeros_like(lon)
            NDFD = pyproj.Proj("+proj=lcc +lat_1=25 +lat_2=25 +lon_0=-95 +x_0=0 +y_0=0 +"
                               "ellps=WGS84 +datum=WGS84 +units=m +no_defs")

            # Convert lat/lon grid into projection space
            X, Y = NDFD(lon, lat)
            WRFX, WRFY = NDFD(wrflon, wrflat)
            # print(WRFX.shape, WRFY.shape)
            # Convert lat/lon reports into projection space
            x, y = NDFD(lons, lats)

            # Create KD-Tree for effecient lookup
            gpoints = np.array(list(zip(X.ravel(), Y.ravel())))
            gtree = sp.spatial.cKDTree(gpoints)

            # Run the KD-Tree to get distances to nearest gridbox and
            # then index of nearest grid point
            dists, inds = gtree.query(np.array(list(zip(x, y))),
                                      distance_upper_bound=1000000000.)

            # Convert index of 1D array into index of 2D lat/lon array
            xind, yind = np.unravel_index(inds[mask], X.shape)

            # Loop through all points and increment that grid cell by 1
            for xi, yi in zip(xind, yind):
                grid[xi, yi] = 1
            # print(grid)

            # Gaussian smoother over our grid to create practically perfect probs
            tmppperf = ndimage.gaussian_filter(grid,sigma=sigma, order=0)

            # Interpolate to WRF grid
            pperf = bilinear_interp(X, Y, WRFX, WRFY, tmppperf)
            print('Min/Max Prac Perf: ', np.min(pperf), '/', np.max(pperf))
        except:
            print("No reports to grid")
            pperf = np.zeros_like(wrflon)

    # Remove report CSV file
    os.remove(rptfile)
    print("Practically Perfect min/max: ", np.min(pperf), ' , ', np.max(pperf))

    return pperf, wrflon, wrflat

def gen_surrogate_severe_reports(uh_arr, sim_refl_arr, uh_thresh,
                                lats, lons, dx=4., spc_grid=False):
    """
    Generates surrogate severe reports given a UH probability array
    and a simulated reflectivity array from the same forecast.

    Inputs
    ------
    uh_arr ------------ 3-D array-like of deterministic updraft helicity
                        forecast values which is structured
                        like so: uh(time, lat, lon)
    sim_refl_arr ------ 4-D array-like of simulated reflectivity to
                        ensure convection is collocated with UH
                        which is structured like so:
                        dbz(time, zlevel, lat, lon)
    uh_thresh --------- float to use as UH threshold for generating
                        SSR's
    lats -------------- 2-D latitude array-like that
                        corresponds to uh_arr and sim_refl_arr
                        grid
    lons -------------- 2-D longitude array-like that
                        corresponds to uh_arr and sim_refl_arr
                        grid
    dx ---------------- horizontal grid-spacing of original
                        dataset as float in kilometers to
                        determine reflectivity mask
    spc_grid ---------- boolean specifying whether to re-grid SSRs
                        to SPC 81-km grid to account for spatial
                        disparity of typical reports

    Outputs
    -------
    returns a 2-D binary array (containing one's in the presence of
    an SSR and zero's otherwise)
    """
    # Initialize surrogate severe report array
    ssr_arr = np.zeros_like(uh_arr[0], dtype=bool)

    # For reflectivity mask
    xinds, yinds = np.meshgrid(np.arange(len(lons[0,:])), np.arange(len(lats[:])))
    # print(np.shape(yinds), np.shape(xinds))

    # Get lats and lons for practically perfect grid
    ppfile = '{}/pperf_grid_template.npz'.format(package_dir)
    f = np.load(ppfile)
    lon = f["lon"]
    lat = f["lat"]
    f.close()

    for ind, hourly_uh_field in enumerate(uh_arr):
        uhmask = (hourly_uh_field >= uh_thresh) # Mask out everything under thresh
        composite_refl = np.max(sim_refl_arr[ind], axis=0) # Find max in column
        uh_inds = np.where(uhmask == True) # Identify UH points
        for fcstind in list(zip(uh_inds[0].ravel(), uh_inds[1].ravel())):
            # Check each UH point for reflectivity > 35.0 dBZ
            r = 25./dx # See Sobash et al 2011; 25-km radius and 35 dBZ
            distmask = dist_mask(fcstind[1], fcstind[0], xinds, yinds, r)
            if (np.asarray(composite_refl)[distmask] < 35.0).all():
                # If no reflectivity within 25 km, negate UH point
                uhmask[fcstind] = False
        # SSR either stays from previous time step, or new SSR from UH/REFL
        ssr_arr = ssr_arr | uhmask

    # If accounting for spatial disparity of storm reports (more realistic),
    #  grid reports to 81-km SPC grid
    if spc_grid:
        # Set up empty grid onto correct projection
        grid = np.zeros_like(lon)
        NDFD = pyproj.Proj("+proj=lcc +lat_1=25 +lat_2=25 +lon_0=-95 +x_0=0 +y_0=0 +"
                           "ellps=WGS84 +datum=WGS84 +units=m +no_defs")

        # Identify SSR locations
        ssr_loc_inds = np.where(ssr_arr == 1)
        # print("Shape of SSR locations", np.shape(ssr_loc_inds))
        # print("SSR location indices", ssr_loc_inds)
        ssr_lons = lons[ssr_loc_inds]
        ssr_lats = lats[ssr_loc_inds]
        # print("SSR Lat/Lons", ssr_lats, ssr_lons)

        # Convert lat/lon grid into projection space
        X, Y = NDFD(lon, lat)
        WRFX, WRFY = NDFD(lons, lats)

        # Convert lat/lon reports into projection space
        x, y = NDFD(ssr_lons, ssr_lats)

        # Create KD-Tree for effecient lookup
        gpoints = np.array(list(zip(X.ravel(), Y.ravel())))
        gtree = sp.spatial.cKDTree(gpoints)

        # Run the KD-Tree to get distances to nearest gridbox and
        # then index of nearest grid point
        dists, inds = gtree.query(np.array(list(zip(x, y))),
                                  distance_upper_bound=1000000000.)

        # Convert index of 1D array into index of 2D lat/lon array
        xind, yind = np.unravel_index(inds, X.shape)

        # Loop through all points and increment that grid cell by 1
        for xi, yi in zip(xind, yind):
            grid[xi, yi] = 1

        # Overwrite SSR array as new 81-km gridded SSRs
        ssr_arr = grid

    return np.array(ssr_arr, dtype=int)

def gen_SSPFs_from_SSRs(ssr_arr, sigma):
    """
    Generates surrogate severe probabilistic forecast from
    gridded surrogate severe reports.

    Inputs
    ------
    ssr_arr ------- 2-D array-like containing binary
                    "surrogate severe reports", where
                    a 1 indicates existance of an SSR
                    and a 0 indicates no SSR.
    sigma --------- float to determine the
                    value of the gaussian smoothing
                    parameter applied to SSRs,
                    equivalent to:
                    neighborhood / horiz-grid spacing

    Outputs
    -------
    returns 2-D array-like of surrogate severe probabilities
    based on the given SSR grid
    """
    return ndimage.gaussian_filter(np.array(ssr_arr, dtype=float),
                                    sigma=sigma, order=0)

def dist_mask(xind, yind, xpts, ypts, r):
    """
    Calculates a mask that evaluates to true in locations
    where gridpoints are within a given radius and false
    in locations where gridpoints are outside the radius.

    Inputs
    ------
    xind ------- the x-index of a gridpoint of interest
    yind ------- the y-index of a gridpoint of interest
    xpts ------- 2D mesh array of x-indices of which to evaluate distances
    ypts ------- 2D mesh array of y-indices of which to evaluate distances
    r ---------- desired radius (in number of grid pts)

    Outputs
    -------
    returns a mask of shape xpts.shape that contains True's
    where the grid is less than the given radius
    """
    return (np.sqrt(((xind - xpts)**2) + ((yind - ypts)**2)) <= r)

def ens_frequency(ens_field, thresh, field='counts', axis=0):
    """
    Calculate ensemble frequency of a given 2-dimensional variable for a specified threshold.

    The number of ensemble members is inferred by the first dimension of `ens_field`,
    unless the axis argument is specified.

    Inputs
    ------
    ens_field --------- P x N x M ndarray; the ensemble field to calculate relative
                        frequencies over. P is the number of ensemble members, N is
                        the number of grid points in the y-dimension, and M is the
                        number of grid points in the x-dimension
    thresh ------------ float; threshold value (greater than or equal to) for which
                        to calculate ensemble relative frequencies.
    field ------------- str; return a numpy ndarray according to the following options:
                        - 'counts' (default): return an N x M array of the ensemble frequency
                        - 'bins': return a P x N x M array of the binary hits and misses for each ensemble member
                        - 'probs': return an N x M array of the ensemble relative frequency
    axis -------------- int; axis of ensemble member dimension. Default is 0

    Outputs
    -------
    (P) x N x M ndarray; the ensemble frequency/binary/probability field(s)
    """
    mask = (ens_field >= thresh)
    members, ypoints, xpoints = np.asarray(mask).nonzero()

    bins = np.zeros_like(ens_field)
    bins[members, ypoints, xpoints] = 1.

    if field == 'bins':
        return bins
    elif field == 'counts':
        return np.sum(bins, axis=axis)
    elif field == 'probs':
        return np.sum(bins, axis=axis) / ens_field.shape[axis]

#############################################################
# Begin verification metrics
#############################################################

def FSS(fcstprobarray, obarray):
    """
    Calculates fractions skill score for a probabilstic
    ensemble forecast, Obs need to be pre-interpolated
    to native model grid. By default, ignores masked
    values

    Inputs
    ------
    fcstprobarray ----- 2-D array-like that contains forecast
                        to be verified.
    obarray ----------- 2-D array-like that contains observations
                        to verify fcstprobarray with. Should be same
                        shape as fcstprobarray.

    Outputs
    -------
    Returns fss as a float for cases in which there are nonzero probabilities for
    one or both the observations and forecast. Otherwise, np.nan is returned.
    """
    # First calculate FBS (Fractions Brier Score) on whole grid
    probs = fcstprobarray
    obs = obarray
    npts = len(fcstprobarray[:,0]) * len(fcstprobarray[0,:])
    fbs = 0.
    fbs_worst = 0.
    print('Max ens probs and max ob probs: ', np.max(probs), np.max(obs))

    # Calculate FBS at each grid point and aggregate.
    if (np.max(obs) > 0.) or (np.max(probs) > 0.):
        for i in range(len(fcstprobarray[0,:])):
            for j in range(len(fcstprobarray[:,0])):
                # print(probs[j,i], obs[j,i])
                if (np.isnan(probs[j,i]) == False) and (np.isnan(obs[j,i]) == False):
                    # print(probs[j,i], obs[j,i])
                    fbs += (probs[j,i] - obs[j,i])**2
                    fbs_worst += probs[j,i]**2 + obs[j,i]**2
        print('FBS: ', fbs)
        fbs, fbs_worst = fbs/npts, fbs_worst/npts
        # Use FBS and FBS worst to calculate FSS for whole grid
        fss = 1 - (fbs/fbs_worst)
        print("FSS and num points: ", fss, ',', npts)
    else:
        print('NULL time/case, cannot calculate FSS')
        return np.nan

    return fss

def FSSnetcdf(probpath, obspath, fhr, var='updraft_helicity',
        thresh=25., rboxpath=None):
    """
    Calculates fractional skill score for a probabilstic
    ensemble forecast from probability and observation
    files. Obs need to be pre-interpolated
    to native model grid.

    Inputs
    ------
    probpath ---------- path to netCDF file containing
                        probability values. IMPORTANT
                        NOTE - prob file is expected to
                        be organized like in probcalcSUBSET.f
                            Variable P_HYD[0,i,:,:]:
                                i         Var
                            [0]   Refl > 40 dBZ probs
                            [1]   UH > 25 m2/s2 probs
                            [2]   UH > 40 m2/s2 probs
                            [3]   UH > 100 m2/s2 probs
                            [4]   Wind Speed > 40 mph probs
    obspath ----------- path to netCDF file containing
                        observation verification values.
                        Currently, practically perfect
                        for verifying UH is the only
                        verification type supported.
    var --------------- string describing variable to verify.
                        Only supports 'updraft_helicity'
                        option as of right now.
    thresh ------------ float describing threshold of
                        variable to use when
                        pulling probs. Choices
                        for UH are 25, 40, and 100 m2/s2.
                        Choices for Reflectivity and
                        Win Speed are 40 (dbz) and
                        40 (mph) respectively.
    rboxpath ---------- optional path to sensitivity calc
                        input file. Only needed if using
                        FSS to verify subsets. Otherwise
                        leave as None.
    prob_var ---------- name of netCDF variable in which
                        probabilities are stored (as string).
                        If using probability calculations from
                        this library, probabilities will by default
                        be stored in 'P_HYD'.

    Outputs
    -------
    returns fss_all, fss_rbox, and sigma back as a
    tuple.

    fss_all -- fractional skill score as float for
                whole domain.
    fss_rbox - fractional skill score as float valid
                for response box. If not verifying
                subsets, then returns 9e9.
    sigma ---- sigma value of practically perfect stored
                in obspath.
    """
    # Open probabilistic forecast and observational datasets
    probdat = Dataset(probpath)
    obsdat = Dataset(obspath)
    initstr = obsdat.START_DATE
    fhrs = obsdat.variables['fhr'][:]
    obind = np.where(fhrs-date2num(fhr, 'hours since ' + initstr) == 0)[0][0]
    print("Fcst hr for FSS calc: ", fhrs[obind])

    # Choose correct indices based on variable and threshold
    probinds = {'reflectivity': {40: 0, 50: 5},
                'updraft_helicity': {25: 1, 40: 2, 100: 3},
                'wind_speed': {40: 4}}
    # If UH, pull practically perfect
    if var == 'updraft_helicity':
        obs = obsdat.variables['practically_perfect'][:]
        sig = obsdat.variables['sigma'][:]
    elif var == 'reflectivity':
        ##### UNDER CONSTRUCTION ########
        pass
    else:
        raise ValueError('Support for {} not yet built in.'.format(var))

    # Pull and splice probability variable
    probvar = probdat.variables['P_HYD'][0]
    d = probinds[var]
    probs = probvar[d[int(thresh)]]

    # If sigma was passed, use it to smooth probs
    # if smooth_w_sigma is not None:
    #     probs = ndimage.gaussian_filter(probs, sigma=smooth_w_sigma)

    # First calculate FBS (Fractions Brier Score) on whole grid
    lats = probdat.variables['XLAT'][0]
    lons = probdat.variables['XLONG'][0]
    npts = len(lats[:,0]) * len(lons[0,:])
    fbs = 0.
    fbs_worst = 0.
    fss_all = 0.
    print('Max ens probs and max ob probs: ', np.max(probs), np.max(obs[obind]))

    # Calculate FBS at each grid point and aggregate.
    if (np.max(obs[obind]) > 0.) or (np.max(probs) > 0.):
        for i in range(len(lons[0,:])):
            for j in range(len(lats[:,0])):
                fbs += (probs[j,i] - obs[obind,j,i])**2
                fbs_worst += probs[j,i]**2 + obs[obind,j,i]**2
        print('FBS: ', fbs)
        fbs, fbs_worst = fbs/npts, fbs_worst/npts
        # Use FBS and FBS worst to calculate FSS for whole grid
        fss_all = 1 - (fbs/fbs_worst)
        print("FSS Total and num points: ", fss_all, ',', npts)

        # Calculate FSS within response box if using ensemble sensitivity
        if rboxpath is not None:
            sensin = np.genfromtxt(rboxpath, dtype=str)
            rbox = sensin[4:8]
            llon, ulon, llat, ulat = np.array(rbox, dtype=float)
            lonmask = (lons > llon) & (lons < ulon)
            latmask = (lats > llat) & (lats < ulat)
            mask = lonmask & latmask
            print(np.shape(mask))
            print("Min/Max Lon:", np.min(lons[mask]), np.max(lons[mask]))
            print("Min/Max Lat:", np.min(lats[mask]), np.max(lats[mask]))
            masked_probs = probs[mask]
            masked_obs = obs[obind][mask]
            npts = len(masked_probs); fbs = 0.; fbs_worst = 0.
            for i in range(len(masked_probs)):
                fbs += (masked_probs[i] - masked_obs[i])**2
                fbs_worst += masked_probs[i]**2 + masked_obs[i]**2
            fbs, fbs_worst = fbs/npts, fbs_worst/npts
            print('FBS: ', fbs)
            print('FBS Worst: ', fbs_worst)
            # Use FBS and FBS worst to calculate FSS for whole grid
            fss_rbox = 1 - (fbs/fbs_worst)
            print("FSS Rbox and num points: ", fss_rbox, ',', npts)
        else:
            fss_rbox = 9e9
    else:
        print('NULL time/case, cannot calculate FSS')

    return fss_all, fss_rbox, sig

def ReliabilityTotal(probpath, runinitdate, fhr, obpath=None, var='updraft_helicity',
                thresh=25., sixhr=False, nbrhd=0.):
    """
    Calculates reliability for a probabilstic
    ensemble forecast over the entire model domain. Obs
    need to be pre-interpolated to native model grid.

    Inputs
    ------
    probpath --- path to netCDF file containing
                 probability values. IMPORTANT
                 NOTE - prob file is expected to
                 be organized like in probcalcSUBSET.f
                  Variable P_HYD[0,i,:,:]:
                   i         Var
                  [0]   Refl > 40 dBZ probs
                  [1]   UH > 25 m2/s2 probs
                  [2]   UH > 40 m2/s2 probs
                  [3]   UH > 100 m2/s2 probs
                  [4]   Wind Speed > 40 mph probs
    obpath ------ path to binary observations gridded to
                  the native WRF domain. If set to None,
                  will automatically calculate gridded obs.
    runinitdate - datetime obj describing model initiation
                  time. Used to pull correct storm reports.
    fhr --------- integer describing response time in number
                  of forecast hours.
    var --------- string describing variable to verify.
                  Only supports 'updraft_helicity'
                  option as of right now.
    thresh ------ float describing threshold of
                  variable to use when
                  pulling probs. Choices
                  for UH are 25, 40, and 100 m2/s2.
                  Choices for Reflectivity and
                  Wind Speed are 40 (dbz) and
                  40 (mph) respectively.
    sixhr ------- option of calculating 1 hours worth
                  or 6 hours worth of SPC storm reports.
                  If set to True, assuming probabilities are
                  also valid over a 6-hour period.
    nbrhd ------- neighborhood distance in km used with the
                  probability calculations.

    Outputs
    -------
    returns an array of probability bins, forecast frequencies for the total domain
    of each bin, and observation hit rates for each bin.
    """
    print('Starting reliability calculations...')
    prob_bins = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # Open probabilistic forecast and observational datasets
    probdat = Dataset(probpath)

    # Choose correct indices based on variable and threshold
    probinds = {'reflectivity' : {40 : 0, 50 : 5},
                'updraft_helicity' : {25 : 1, 40 : 2, 100 : 3},
                'wind_speed' : {40 : 4}}

    if var == 'updraft_helicity':
        if obpath is not None:
            dat = Dataset(obpath)
            times = dat.variables['fhr'][:]
            inds = np.where(times == fhr)
            grid = dat.variables['nearest_neighbor'][inds][0]
        else:
            grid = nearest_neighbor_ttu(runinitdate, sixhr, fhr, nbrhd=0.)
    else:
        raise ValueError('Sorry, support for {} is not yet built in.'.format(var))
    # Pull and splice probability variable
    probvar = probdat.variables['P_HYD'][0]
    d = probinds[var]
    fcstprobs = probvar[d[int(thresh)]]
    wrf.disable_xarray()
    lats = wrf.getvar(probdat, 'lat')
    lons = wrf.getvar(probdat, 'lon')
    dx = probdat.DX / 1000. # dx in km
    # Get arrays of x and y indices for distance calculations
    xinds, yinds = np.meshgrid(np.arange(len(lons[0,:])), np.arange(len(lats[:])))

    # Sort probabilities into bins
    fcstfreq_tot = np.zeros((len(prob_bins)))   # N probs falling into bin for whole domain
    ob_hr_tot = np.zeros((len(prob_bins)))      # Ob hit rate for bin and whole domain

    for i in range(len(prob_bins)):
        hits = 0
        prob = prob_bins[i]
        fcstinds = np.where((np.abs(fcstprobs - prob) <= 10) & (fcstprobs < prob))
        fcstfreq_tot[i] = len(fcstinds[0])
        if fcstfreq_tot[i] == 0:
            ob_hr_tot[i] = 0
        else:
            for fcstind in list(zip(fcstinds[0].ravel(), fcstinds[1].ravel())):
                r = nbrhd/dx
                mask = dist_mask(fcstind[1], fcstind[0], xinds, yinds, r)
                masked_grid = np.ma.masked_array(grid, mask=~mask)
                if (masked_grid == 1).any():
                    hits += 1
            tot =  len(list(zip(fcstinds[0].ravel(),
                                        fcstinds[1].ravel())))
            ob_hr_tot[i] = hits / tot
            print("Total Hits/Tot: ", hits, tot)

    return prob_bins, fcstfreq_tot, ob_hr_tot

def ReliabilityRbox(probpath, runinitdate, fhr,  rboxpath, obpath=None,
                var='updraft_helicity',
                thresh=25.,sixhr=False, nbrhd=0.):
    """
    Calculates reliability for a probabilstic
    ensemble forecast and returns it. Obs
    need to be pre-interpolated to native model grid.

    Inputs
    ------
    probpath --- path to netCDF file containing
                 probability values. IMPORTANT
                 NOTE - prob file is expected to
                 be organized like in probcalcSUBSET.f
                  Variable P_HYD[0,i,:,:]:
                   i         Var
                  [0]   Refl > 40 dBZ probs
                  [1]   UH > 25 m2/s2 probs
                  [2]   UH > 40 m2/s2 probs
                  [3]   UH > 100 m2/s2 probs
                  [4]   Wind Speed > 40 mph probs
    obpath ------ path to binary observations gridded to
                  the native WRF domain. If set to None,
                  will automatically calculate gridded obs.
    runinitdate - datetime obj describing model initiation
                  time. Used to pull correct storm reports.
    fhr --------- integer describing response time in number
                  of forecast hours.
    var --------- string describing variable to verify.
                  Only supports 'updraft_helicity'
                  option as of right now.
    thresh ------ float describing threshold of
                  variable to use when
                  pulling probs. Choices
                  for UH are 25, 40, and 100 m2/s2.
                  Choices for Reflectivity and
                  Wind Speed are 40 (dbz) and
                  40 (mph) respectively.
    rboxpath ---- absolute path to sensitivity calc
                  input file. Only needed if verifying
                  subsets. Otherwise
                  leave as None.
    sixhr ------- option of calculating 1 hours worth
                  or 6 hours worth of SPC storm reports.
                  If set to True, assuming probabilities are
                  also valid over a 6-hour period.
    nbrhd ------- neighborhood distance in km used with the
                  probability calculations.

    Outputs
    -------
    returns an array of probability bins, forecast frequencies for the total domain
    of each bin, observation hit rates for each bin, forecast frequencies for the
    response box for each bin, observation hit rates for the response box for each
    bin (if not verifying subsets, will return arrays of zeros for the response box
    metrics).
    """
    print('Starting reliability calculations...')
    prob_bins = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # Open probabilistic forecast and observational datasets
    probdat = Dataset(probpath)

    # Choose correct indices based on variable and threshold
    probinds = {'reflectivity' : {40 : 0, 50 : 5},
                'updraft_helicity' : {25 : 1, 40 : 2, 100 : 3},
                'wind_speed' : {40 : 4}}

    if var == 'updraft_helicity':
        if obpath is not None:
            dat = Dataset(obpath)
            times = dat.variables['fhr'][:]
            inds = np.where(times == fhr)
            grid = dat.variables['nearest_neighbor'][inds][0]
        else:
            grid = nearest_neighbor_ttu(runinitdate, sixhr, fhr, nbrhd=0.)
    else:
        raise ValueError('Sorry, support for {} is not yet built in.'.format(var))
    # Pull and splice probability variable
    probvar = probdat.variables['P_HYD'][0]
    d = probinds[var]
    fcstprobs = probvar[d[int(thresh)]]
    wrf.disable_xarray()
    lats = wrf.getvar(probdat, 'lat')
    lons = wrf.getvar(probdat, 'lon')
    dx = probdat.DX / 1000. # dx in km
    r = nbrhd/dx
    # Get arrays of x and y indices for distance calculations
    xinds, yinds = np.meshgrid(np.arange(len(lons[0,:])), np.arange(len(lats[:])))

    # Sort probabilities into bins
    fcstfreq_tot = np.zeros((len(prob_bins)))   # N probs falling into bin for whole domain
    fcstfreq_rbox = np.zeros((len(prob_bins)))  # N probs in rbox falling into bin
    ob_hr_tot = np.zeros((len(prob_bins)))      # Ob hit rate for bin and whole domain
    ob_hr_rbox = np.zeros((len(prob_bins)))     # Ob hit rate for bin in rbox

    # Create mask for isolating response box in grid if applicable
    if rboxpath is not None:
        sensin = np.genfromtxt(rboxpath, dtype=str)
        rbox = sensin[4:8]
        llon, ulon, llat, ulat = np.array(rbox, dtype=float)
        lonmask = (lons > llon) & (lons < ulon)
        latmask = (lats > llat) & (lats < ulat)
        mask = lonmask & latmask
        # print(len(lons[xinds[mask]]))
        print(fcstprobs.shape)
        masked_probs = np.ma.masked_array(fcstprobs, mask=~mask)
        print(llat, ulat, llon, ulon)
        print(lats[np.ma.where(masked_probs)])
        print(np.ma.where(masked_probs))
        print("Number of Rbox points:", len(np.ma.where(masked_probs)[0]))

    for i in range(len(prob_bins)):
        hits = 0
        prob = prob_bins[i]
        print("Calculating reliability over rbox for prob bin {}-{}%".format(prob-10,
                                                                    prob))
        # If verifying subsets, we want the reliability inside the response box
        if rboxpath is not None:
            fcstinds = np.ma.where((np.abs(masked_probs - prob) <= 10) & \
                                (masked_probs < prob))
            # print(fcstinds)
            fcstfreq_rbox[i] = len(fcstinds[0])
            if fcstfreq_rbox[i] == 0:
                ob_hr_rbox[i] = 0
            else:
                import matplotlib.pyplot as plt
                # fig, ((ax)) = plt.subplots(1, 1, figsize=(15,15))
                for fcstind in list(zip(fcstinds[0].ravel(), fcstinds[1].ravel())):
                    # print(lats[fcstind[0], fcstind[1]], lons[fcstind[0], fcstind[1]])
                    mask = dist_mask(fcstind[0], fcstind[1], yinds, xinds, r)
                    masked_grid = np.ma.masked_array(grid, mask=~mask)
                    # Dost thou desire a sanity check?
                    # clevs = np.linspace(0,2,4)
                    # ax.contourf(grid, clevs)
                    # ax.contourf(masked_grid, clevs, cmap="jet")
                    if (masked_grid == 1).any():
                        hits += 1
                tot = len(list(zip(fcstinds[0].ravel(),
                                            fcstinds[1].ravel())))
                ob_hr_rbox[i] = hits / tot
                print("Rbox Hits/Total: ", hits, tot)
                # ax.set_title("Fcst Frequency: {}; Ob Hit Rate: {}".format(fcstfreq_rbox[i],
                                    # ob_hr_rbox[i]))
                # plt.savefig("slow_rel_prob{}-{}.png".format(prob-10,prob))

    return prob_bins, fcstfreq_rbox, ob_hr_rbox

def rmse(predictions, targets, axis=None, nan=False):
    """
    Root Mean Square Error (RMSE)

    Calculate RMSE on grid or timeseries

    Inputs
    ------
    predictions --- array, forecast variable array.
    targets ------- array, observation variable array.
    axis ---------- tuple, optional. Axes over which to perform calculation, If none, RMSE
                    calculated for entire array shape.
    nan ----------- bool, optional. Determines if nan values in inputs should result in
                    nan output. Default is to ignore nan points.
    Outputs
    -------
    Root Mean Square Error of predictions
    """
    if nan:
        rmse_data = np.sqrt(np.mean(((predictions - targets) ** 2), axis=axis))
    else:
        rmse_data = np.sqrt(np.nanmean(((predictions - targets) ** 2), axis=axis))
    return rmse_data

def mse(predictions, targets, axis=None, nan=False):
    """
    Mean Square Error (MSE)

    Calculate MSE on grid or timeseries

    Inputs
    ------
    predictions --- array, forecast variable array.
    targets ------- array, observation variable array.
    axis ---------- tuple, optional. Axes over which to perform calculation, If none, RMSE
                    calculated for entire array shape.
    nan ----------- bool, optional. Determines if nan values in inputs should result in
                    nan output. Default is to ignore nan points.
    Outputs
    -------
    Mean Square Error of predictions
    """
    if nan:
        mse_data = np.mean(((predictions - targets) ** 2), axis=axis)
    else:
        mse_data = np.nanmean(((predictions - targets) ** 2), axis=axis)
    return mse_data

def mae(predictions, targets, axis=None, nan=False):
    """
    Mean Absolute Error (MAE)

    Calculates MAE on grid or timeseries

    Inputs
    ------
    predictions --- array, forecast variable array.
    targets ------- array, observation variable array.
    axis ---------- tuple, optional. Axes over which to perform calculation.
                    If none, MAE calculated over entire array.
    nan ----------- bool, optional. Determines if nan values in inputs should
                    result in nan output. Default is to ignore nan points.

    Outputs
    -------
    Returns mean absolute error of predictions
    """
    if nan:
        mae_data = np.mean((np.abs(predictions - targets)), axis=axis)
    else:
        mae_data = np.nanmean((np.abs(predictions - targets)), axis=axis)
    return mae_data

######### GridRad Reflectivity Verification Metrics ##################
def calc_subset_avg_response_rbox(member_list, rvalues_ncfile, rfuncstr):
    """
    Calculates the average response function within a
    response box over a specified set of one-based subset members.

    Inputs
    ------
    member_list ----------- list of one-based integers specifying subset
                            members to use in calculation
    rvalues_ncfile -------- absolute path to netCDF file containing
                            response function (reflectivity) values
                            from ensemble members. Should be output
                            of a sixhresens.f call
    rfuncstr -------------- response function string specifying which
                            variable to average. Can be any of the following:
                            (1) DBZ_AVG
                            (2) DBZ_MAX
                            (3) UH_AVG
                            (4) UH_MAX
                            (5) PCP
                            (6) WSPD_AVG
                            (7) UH_COV
                            (8) DBZ_COV

    Outputs
    -------
    returns the average response function value valid in the response
    box over the specified ensemble members.
    """
    ds = xr.open_dataset(rvalues_ncfile)
    rfunc_vals = ds[rfuncstr].values
    return np.nanmean(rfunc_vals[np.asarray(member_list, dtype=int)-1])

def calc_refl_max_rbox(gridradfiles, rboxpath, zlev):
    """
    Calculates the reflectivity maximum value in a given response
    box. This function is meant to be used to verify ESA-based
    subsets that use reflectivity maxima as their response function.

    Inputs
    ------
    gridradfiles ---------- list of absolute paths to GridRad file
    rboxpath -------------- absolute path to sensitivity input file for
                            sensitivity calculations that contains the
                            response box bounds
    zlev ------------------ altitude index to slice reflectivity from

    Outputs
    -------
    returns the response box reflectivity maximum along with its
    lat / lon location as a tuple like so: (refl_max, lat, lon)
    """
    # Read response box bounds
    esensin = np.genfromtxt(rboxpath)
    rbox_bounds = esensin[4:8]

    dat = pp.merge_refl_data(gridradfiles)
    lons = dat["Longitude"]
    lats = dat["Latitude"]

    # Build rbox mask
    lonmask = (lons >= rbox_bounds[0]) & (lons < rbox_bounds[1])
    latmask = (lats >= rbox_bounds[2]) & (lats < rbox_bounds[3])
    mask = latmask & lonmask

    # Find reflectivity maximum and its location
    first_refl = dat.values[0,zlev]
    refl_max = np.nanmax(first_refl[mask])
    for i in range(len(dat.values)):
        refl_tmp = np.nanmax(dat.values[i,zlev][mask])
        if refl_tmp > refl_max:
            refl_max = refl_tmp
    # reflmaxind = np.where(dat.values == refl_max)
    # meshlats, meshlons = np.meshgrid(lats, lons)
    # reflmax_lon = meshlons[reflmaxind[2:]][0]
    # reflmax_lat = meshlats[reflmaxind[2:]][0]

    return refl_max

def calc_refl_cov_rbox(interpgridradfiles, rboxpath, zlev, refl_thresh):
    """
    Calculates the reflectivity coverage (in number of grid points)
    over a given response box. This function is meant to be used to
    verify ESA-based subsets that use reflectivity coverage as their
    response function.

    Inputs
    ------
    interpgridradfiles ---- absolute path to original GridRad file
    rboxpath -------------- absolute path to sensitivity input file for
                            sensitivity calculations that contains the
                            response box bounds
    zlev ------------------ zero-based integer describing which altitude
                            level to pull reflectivity data from
    refl_thresh ----------- float describing reflectivity threshold to
                            calculate coverage with

    Outputs
    -------
    returns the response box reflectivity coverage as an integer
    as well as the number of grid points in the response box
    """
    # Read response box bounds
    esensin = np.genfromtxt(rboxpath)
    rbox_bounds = esensin[4:8]

    # Read data
    dats = [xr.open_dataset(file) for file in interpgridradfiles]
    print(dats[0])
    time = [dats[0].STARTDATE]
    lats = dats[0]["XLAT"][0]
    lons = dats[0]["XLONG"][0]
    refl = dats[0]["GridRad_Refl"].values
    dat = xr.DataArray(refl)#, coords=[time, lats, lons],
                        #dims=['times', 'latitude', 'longitude'])
    dat.name = "Reflectivity"
    for i in range(1,len(dats)):
        # Pull reflectivity information
        ds = dats[i]
        time = [ds.STARTDATE]
        lats = ds["XLAT"][0]
        lons = ds["XLONG"][0]
        refl = ds["GridRad_Refl"].values
        da = xr.DataArray(refl)#, coords=[time, lats, lons],
                                #dims=['times', 'latitude', 'longitude'])
        da.name = "Reflectivity"
        dat = xr.concat([dat, da], dim='times')

    # Build rbox mask
    lonmask = (lons >= rbox_bounds[0]) & (lons < rbox_bounds[1])
    latmask = (lats >= rbox_bounds[2]) & (lats < rbox_bounds[3])
    mask = latmask & lonmask

    # Capture grid points where reflectivity exceeds threshold
    refl_cov = 0
    for i in range(len(dat.values)):
        refl_rbox = dat.values[i,zlev][mask]
        refl_exceeds = refl_rbox[refl_rbox > refl_thresh]
        refl_cov += len(refl_exceeds)
    npts_rbox = len(refl_rbox)

    return refl_cov, npts_rbox

def calc_og_refl_cov_rbox(gridradfiles, rboxpath, zlev, refl_thresh):
    """
    Calculates the reflectivity coverage (in number of grid points)
    over a given response box. This function is meant to be used to
    verify ESA-based subsets that use reflectivity coverage as their
    response function.

    Inputs
    ------
    gridradfiles ---------- absolute path to original GridRad file
    rboxpath -------------- absolute path to sensitivity input file for
                            sensitivity calculations that contains the
                            response box bounds
    zlev ------------------ zero-based integer describing which altitude
                            level to pull reflectivity data from
    refl_thresh ----------- float describing reflectivity threshold to
                            calculate coverage with

    Outputs
    -------
    returns the response box reflectivity coverage as an integer
    as well as the number of grid points in the response box
    """
    # Read response box bounds
    esensin = np.genfromtxt(rboxpath)
    rbox_bounds = esensin[4:8]

    # Read data
    dats = [xr.open_dataset(file) for file in gridradfiles]
    inds = dats[0]["index"]
    time = dats[0].Analysis_time
    lats = dats[0]["Latitude"]
    lons = dats[0]["Longitude"]-360.
    z = dats[0]["Altitude"]
    refl = dats[0]["Reflectivity"].values
    refl_vals = np.zeros(len(z.values)*len(lats.values)*len(lons.values))*np.NaN
    refl_vals[inds.values] = refl
    refl_reshape = refl_vals.reshape((len(z.values), len(lats.values),
                                        len(lons.values)))
    dat = xr.DataArray(refl_reshape, coords=[z, lats, lons])
    dat.name = "Reflectivity"
    for i in range(1,len(dats)):
        # Pull reflectivity information
        ds = dats[i]
        inds = ds["index"]
        time = ds.Analysis_time
        lats = ds["Latitude"]
        lons = ds["Longitude"]-360.
        z = ds["Altitude"]
        refl = ds["Reflectivity"].values
        refl_vals = np.zeros(len(z.values)*len(lats.values)*len(lons.values))*np.NaN
        refl_vals[inds.values] = refl
        refl_reshape = refl_vals.reshape((len(z.values), len(lats.values),
                                            len(lons.values)))
        da = xr.DataArray(refl_reshape, coords=[z, lats, lons])
        da.attrs['Analysis_Endtime'] = time
        da.name = "Reflectivity"
        dat = xr.concat([dat, da], dim='Hours')

    # Build rbox mask
    lonmask = (lons >= rbox_bounds[0]) & (lons < rbox_bounds[1])
    latmask = (lats >= rbox_bounds[2]) & (lats < rbox_bounds[3])
    mask = latmask & lonmask

    # Capture grid points where reflectivity exceeds threshold
    refl_cov = 0
    for i in range(len(dat.values)):
        refl_rbox = dat.values[i,zlev][mask]
        refl_exceeds = refl_rbox[refl_rbox > refl_thresh]
        refl_cov += len(refl_exceeds)
    npts_rbox = len(refl_rbox)

    return refl_cov, npts_rbox
