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
from netCDF4 import Dataset, date2num, num2date
from datetime import datetime, timedelta
from subprocess import call
from scipy import ndimage
import matplotlib.patches as patches
import cartopy.crs as ccrs
import pyproj
import scipy as sp
import os
import csv
import sys

def bilinear_interp(grid1x, grid1y, grid2x, grid2y, z):
    '''
    A method which interpolates a function 
    z(grid1x, grid1y) of a grid (grid1x, grid1y) to another 
    grid (grid2x, grid2y). Returns an array from the approximated
    function of the second grid (approximation of z(grid2x, grid2y)).
    '''
    # Pair flattened x and y values as coordinates
    coords_from = list(zip(grid1y.flatten(), grid1x.flatten()))
    Z = z.flatten()
    # Set up interpolation function with original grid and interp variable
    interp = sp.interpolate.LinearNDInterpolator(coords_from, Z, fill_value=9e9)
    # Interpolate to new grid
    interpolated_z = interp(grid2y, grid2x)
    
    return interpolated_z

def calc_prac_perf(runinitdate, sixhr, rtime, sigma=2):
    '''
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
    sigma -------- optional integer specifying sigma
                    to use for Gaussian filter.
                    Operational practically perfect uses
                    the default sgma of 2.
    Outputs
    -------
    returns tuple containing pract perf probs 
    and lons/lats (respectively) of pperf grid
    '''
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
    
    #Get yesterday's reports CSV file from web
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
    	
    #Get lats and lons for practically perfect grid
    ppfile = '/lustre/work/aucolema/scripts/pperf_grid_template.npz'
    f = np.load(ppfile)
    lon = f["lon"]
    lat = f["lat"]
    f.close()
    
    # Get WRF lats/lons as pperf grid
    ppfile = '/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREFd2'
    dat = Dataset(ppfile)
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
        else:
            hour = rdate.hour
            mask = [(hr == (hour-1)%24) for hr in time]
            print('Reports valid {} to {}'.format((hour-1)%24,hour))
            #print('Report hours from reports used:', np.array(time)[mask])
            #print(mask)
            #print(time)
        try:
            #Set up empty grid onto correct projection
            grid = np.zeros_like(lon)
            NDFD = pyproj.Proj("+proj=lcc +lat_1=25 +lat_2=25 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
        
            # Convert lat/lon grid into projection space
            X, Y = NDFD(lon, lat) 
            WRFX, WRFY = NDFD(wrflon, wrflat)
            #print(WRFX.shape, WRFY.shape)
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
            #print(grid)
    
        	# Gaussian smoother over our grid to create practically perfect probs
            tmppperf = ndimage.gaussian_filter(grid,sigma=sigma, order=0)
            
            # Interpolate to WRF grid
            pperf = bilinear_interp(X, Y, WRFX, WRFY, tmppperf)
            #print(np.shape(pperf))
        except:
            pperf = np.zeros_like(wrflon)
    
    #Remove report CSV file
    os.remove(rptfile)    
    #print(sigma)
    print("Practically Perfect min/max: ", np.min(pperf), ' , ', np.max(pperf))
    
    return pperf, wrflon, wrflat
    
#############################################################
# Begin verification metrics
#############################################################    
    
def FSS(probpath, obspath, fhr, var='updraft_helicity', 
        thresh=25., rboxpath=None):
    '''
    Calculates fractional skill score for a probabilstic
    ensemble forecast and stores it in a csv file. Obs
    need to be pre-interpolated to native model grid.
    
    Inputs
    ------
    probpath - path to netCDF file containing 
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
    obspath -- path to netCDF file containing
                observation verification values.
                Currently, practically perfect
                for verifying UH is the only
                verification type supported.
    var ------ string describing variable to verify.
                Only supports 'updraft_helicity'
                option as of right now.
    thresh --- float describing threshold of
                variable to use when 
                pulling probs. Choices
                for UH are 25, 40, and 100 m2/s2.
                Choices for Reflectivity and
                Win Speed are 40 (dbz) and
                40 (mph) respectively.
    rboxpath - optional path to sensitivity calc
                input file. Only needed if using
                FSS to verify subsets. Otherwise
                leave as None.
    
    Outputs
    -------
    Returns fss_all, fss_rbox, and sigma back as a 
    tuple. 
    
    fss_all -- fractional skill score as float for
                whole domain.
    fss_rbox - fractional skill score as float valid
                for response box. If not verifying 
                subsets, then returns 9e9.
    sigma ---- same sigma specified by input for
                keeping track.
    '''
    # Open probabilistic forecast and observational datasets
    probdat = Dataset(probpath)
    obsdat = Dataset(obspath)
    initstr = obsdat.START_DATE
    fhrs = obsdat.variables['fhr'][:]
    obind = np.where(fhrs-date2num(fhr, 'hours since ' + initstr) == 0)[0][0]
    print("Fcst hr for FSS calc: ", fhrs[obind])

    # Choose correct indices based on variable and threshold
    probinds = {'reflectivity' : {40 : 0}, 
                'updraft_helicity' : {25 : 1, 40 : 2, 100 : 3},
                'wind_speed' : {40 : 4}}
    # If UH, pull practically perfect
    if var == 'updraft_helicity':
        obs = obsdat.variables['practically_perfect'][:]
        sig = obsdat.variables['sigma'][:]
    else:
        raise ValueError('Support for {} not yet built in.'.format(var))
            
    # Pull and splice probability variable
    probvar = probdat.variables['P_HYD'][0]
    d = probinds[var]
    probs = probvar[d[int(thresh)]]
    # First calculate FBS (Fractions Brier Score) on whole grid
    wrf.disable_xarray()
    lats = wrf.getvar(probdat, 'lat')
    lons = wrf.getvar(probdat, 'lon')
    npts = len(lats[:,0]) * len(lons[0,:])
    #prob_gt_than = (probs >= probthresh)
    fbs = 0.
    fbs_worst = 0.
    print('Max ens probs and max ob probs: ', np.max(probs), np.max(obs[obind]))
    if (np.max(obs[obind]) > 0.) or (np.max(probs) > 0.):
        #print(np.shape(probs), np.shape(obs[obind]))
        for i in range(len(lons[0,:])):
            for j in range(len(lats[:,0])):
                fbs += (probs[j,i] - obs[obind,j,i])**2
                fbs_worst += probs[j,i]**2 + obs[obind,j,i]**2
                #print(probs[j,i], obs[obind,j,i])
        print('FBS: ', fbs)
        fbs, fbs_worst = fbs/npts, fbs_worst/npts
        # Use FBS and FBS worst to calculate FSS for whole grid
        fss_all = 1 - (fbs/fbs_worst)
        print("FSS Total and num points: ", fss_all, ',', npts)
        
        # Calculate FSS within response box if using sensitivity
        if rboxpath is not None:
            sensin = np.genfromtxt(rboxpath, dtype=str)
            rbox = sensin[4:8]
            llon, ulon, llat, ulat = np.array(rbox, dtype=float)
            lonmask = (lons > llon) & (lons < ulon)
            latmask = (lats > llat) & (lats < ulat)
            mask = lonmask & latmask
            masked_probs = probs[mask]
            masked_obs = obs[obind][mask]
            npts = len(probs[mask])
            fbs = 0.
            fbs_worst = 0.
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
    
        
