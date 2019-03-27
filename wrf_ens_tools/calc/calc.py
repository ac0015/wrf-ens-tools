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
from profilehooks import profile

package_dir =  os.path.dirname(os.path.abspath(__file__))

@profile
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


def nearest_neighbor_spc(runinitdate, sixhr, rtime, nbrhd=0.,
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


def calc_prac_perf(runinitdate, sixhr, rtime, sigma=2):
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
            for xi, yi in zip(xind, yind):
                grid[xi, yi] = 1

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
            pperf = np.zeros_like(wrflon)

    # Remove report CSV file
    os.remove(rptfile)
    print("Practically Perfect min/max: ", np.min(pperf), ' , ', np.max(pperf))

    return pperf, wrflon, wrflat


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

#############################################################
# Begin verification metrics
#############################################################

def FSS(probpath, obspath, fhr, var='updraft_helicity',
        thresh=25., rboxpath=None):
    """
    Calculates fractional skill score for a probabilstic
    ensemble forecast, Obs need to be pre-interpolated
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
    smooth_w_sigma ---- optional smoothing parameter. If you
                        want to smooth the ensemble probs,
                        replace None default with a sigma
                        value for the standard deviation of
                        the Gaussian kernel to use. Otherwise
                        FSS is calculated with the raw ensemble
                        probabilities.

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
    probinds = {'reflectivity': {40: 0},
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
    wrf.disable_xarray()
    lats = wrf.getvar(probdat, 'lat')
    lons = wrf.getvar(probdat, 'lon')
    npts = len(lats[:,0]) * len(lons[0,:])
    fbs = 0.
    fbs_worst = 0.
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
            print("Min/Max Lon:", np.min(lons[lonmask]), np.max(lons[lonmask]))
            print("Min/Max Lat:", np.min(lats[latmask]), np.max(lats[latmask]))
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
    probinds = {'reflectivity' : {40 : 0},
                'updraft_helicity' : {25 : 1, 40 : 2, 100 : 3},
                'wind_speed' : {40 : 4}}

    if var == 'updraft_helicity':
        if obpath is not None:
            dat = Dataset(obpath)
            times = dat.variables['fhr'][:]
            inds = np.where(times == fhr)
            grid = dat.variables['nearest_neighbor'][inds][0]
        else:
            grid = nearest_neighbor_spc(runinitdate, sixhr, fhr, nbrhd=0.)
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
    yinds, xinds = np.meshgrid(np.arange(len(lons[0,:])), np.arange(len(lats[:])))

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

@profile
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
    probinds = {'reflectivity' : {40 : 0},
                'updraft_helicity' : {25 : 1, 40 : 2, 100 : 3},
                'wind_speed' : {40 : 4}}

    if var == 'updraft_helicity':
        if obpath is not None:
            dat = Dataset(obpath)
            times = dat.variables['fhr'][:]
            inds = np.where(times == fhr)
            grid = dat.variables['nearest_neighbor'][inds][0]
        else:
            grid = nearest_neighbor_spc(runinitdate, sixhr, fhr, nbrhd=0.)
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

@profile
def scipyReliabilityRbox(probpath, runinitdate, fhr,  rboxpath,
                obpath=None, var='updraft_helicity',
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
    probinds = {'reflectivity' : {40 : 0},
                'updraft_helicity' : {25 : 1, 40 : 2, 100 : 3},
                'wind_speed' : {40 : 4}}

    if var == 'updraft_helicity':
        if obpath is not None:
            dat = Dataset(obpath)
            times = dat.variables['fhr'][:]
            inds = np.where(times == fhr)
            grid = dat.variables['nearest_neighbor'][inds][0]
        else:
            grid = nearest_neighbor_spc(runinitdate, sixhr, fhr, nbrhd=0.)
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
    yinds, xinds = np.meshgrid(np.arange(len(lons[0,:])), np.arange(len(lats[:])))

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
        mask = latmask & lonmask
        masked_probs = np.ma.masked_array(fcstprobs, mask=~mask)

    # Apply smoother to act as a neighborhood
    # nbrhd_grid = ndimage.gaussian_filter(grid, sigma=r, order=0)
    # print(grid[grid>0])
    # print(nbrhd_grid[grid>0])
    import matplotlib.pyplot as plt
    # fig, ((ax1, ax2)) = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(15,15))
    # clevs = np.linspace(0,1,100)
    # og = ax1.contourf(grid, clevs, cmap="jet")
    # plt.colorbar(og, ax=ax1)
    # clevs = np.linspace(0,0.5,100)
    # smooth = ax2.contourf(nbrhd_grid, clevs, cmap="jet")
    # plt.colorbar(smooth, ax=ax2)
    # clevs = np.linspace(-0.3,0.1,100)
    # # print(np.min(nbrhd_grid-grid), np.max(nbrhd_grid-grid))
    # # diff = ax3.contourf(nbrhd_grid-grid, clevs, cmap="viridis")
    # # plt.colorbar(diff, ax=ax3)
    # plt.savefig("scipy_smoothed_grid.png")

    for i in range(len(prob_bins)):
        hits = 0
        prob = prob_bins[i]
        print("Calculating reliability over rbox for prob bin {}-{}%".format(prob-10,
                                                                    prob))
        # If verifying subsets, we want the reliability inside the response box
        if rboxpath is not None:
            # Find indices where fcst probs fall into bin
            # fcstinds = np.ma.where((np.abs(masked_probs - prob) <= 10) & \
            #                     (masked_probs < prob))
            fcstmask = (np.abs(fcstprobs - prob) <= 10) & (fcstprobs < prob)
            fcstmask = fcstmask.astype(float)
            print(fcstmask[mask&(fcstmask>0)])
            fuzzy_fcstmask = ndimage.uniform_filter(fcstmask, size=r)
            fuzzy_fcstmask = (fuzzy_fcstmask > 0)
            print(fuzzy_fcstmask[fuzzy_fcstmask == True])
            # print(fuzzy_fcstmask[fuzzy_fcstmask>0.0000000001 & mask])
            combo_mask = (fuzzy_fcstmask > 0)
            masked_grid = np.ma.masked_array(grid, mask=~(mask & fuzzy_fcstmask))
            fig = plt.figure(figsize=(15,15))
            clevs = np.linspace(0,1,50)
            yo = plt.contourf(masked_grid, clevs, cmap='jet')
            plt.colorbar(yo)
            fcstfreq_rbox[i] = len(fuzzy_fcstmask[(mask & fuzzy_fcstmask)])
            if fcstfreq_rbox[i] == 0:
                ob_hr_rbox[i] = 0
            else:
                # Subset of smoothed ob grid that exceeds 0 constitutes a hit
                where = np.ma.where(masked_grid[mask & (fuzzy_fcstmask)] > 0)
                hits = len(where[0])
                # Define ob hit rate
                ob_hr_rbox[i] = hits / fcstfreq_rbox[i]
                print("Rbox Hits/Total: ", hits, fcstfreq_rbox[i])
            plt.title("Fcst Freq: {}; Hits: {}".format(fcstfreq_rbox[i],
                        hits))
            plt.savefig("scipy_prob{}-{}.png".format(prob-10,prob))
            plt.close()

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

######### GridRad Reflectivity Verification Metrics ##################
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
    reflmaxind = np.where(dat.values == refl_max)
    meshlats, meshlons = np.meshgrid(lats, lons)
    reflmax_lon = meshlons[reflmaxind[2:]][0]
    reflmax_lat = meshlats[reflmaxind[2:]][0]

    return refl_max, reflmax_lat, reflmax_lon

def calc_refl_cov_rbox(gridradfiles, rboxpath, zlev, refl_thresh):
    """
    Calculates the reflectivity coverage (in number of grid points)
    over a given response box. This function is meant to be used to
    verify ESA-based subsets that use reflectivity coverage as their
    response function.

    Inputs
    ------
    gridradfiles ---------- absolute path to interpolated GridRad file
    rboxpath -------------- absolute path to sensitivity input file for
                            sensitivity calculations that contains the
                            response box bounds
    zlev ------------------ zero-based integer describing which altitude
                            level to pull reflectivity data from
    refl_thresh ----------- float describing reflectivity threshold to
                            that

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
