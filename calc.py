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

def calc_prac_perf(runinitdate, sixhr, rtime):
    '''
    Implementation of SPC practically perfect
    calculations adapted from SPC code
    
    Inputs
    ------
    runinitdate - datetime obj for 
                    SPC storm reports
    sixhr-------- boolean specifying whether
                    to calculate practically perfect
                    probs over six hr time window or 
                    use one hr time window
    rtime ------- time (in num fcst hrs
                    from runinit) to obtain six hr
                    practically perfect
    Outputs
    -------
    returns tuple containing pract perf probs 
    and lons/lats (respectively) of pperf grid
    '''
    #Get initialization date
    runinitdatef = runinitdate.strftime('%y%m%d')
    
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
    zppfile = 'pperf_grid_template.npz'
    f = np.load(ppfile)
    lon = f["lon"]
    lat = f["lat"]
    print(np.shape(lon))
    f.close()
    
    # Get WRF lats/lons as pperf grid
    #ppfile = '/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREFd2'
    #dat = Dataset(ppfile)
    #lon = dat.variables['XLONG'][0]
    #lat = dat.variables['XLAT'][0]
    #print(np.shape(lon))
    #dat.close()
    
    #If there aren't any reports, practically perfect is zero across grid
    if length == 0:
    	pperf = np.zeros_like(lon)    
    #Otherwise, let's grid the reports
    else:
        # If six hour, mask reports by valid times in window
        if sixhr:
            rdate = runinitdate + timedelta(hours=rtime)
            hour = rdate.hour
            hours = [(hour - i)%24 for i in range(1,7)]
            mask = [(hr in hours) for hr in time]
        else:
            mask = [(hr == (rtime-1)%24) for hr in time]
            #print(mask)
            #print(time)
        try:
            #Set up empty grid onto correct projection
            grid = np.zeros_like(lon)
            NDFD = pyproj.Proj("+proj=lcc +lat_1=25 +lat_2=25 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
        
            # Convert lat/lon grid into projection space
            X, Y = NDFD(lon, lat) 
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
    
        	#Gaussian smoother over our grid to create practically perfect probs
            if sixhr:
                pperf = ndimage.gaussian_filter(grid,sigma=1,order=0)
            else:
                pperf = ndimage.gaussian_filter(grid,sigma=0,order=0)
        except:
            pperf = np.zeros_like(lon)
    
    #Remove report CSV file
    os.remove(rptfile)
    
    return pperf, lon, lat

def FSS(probpath, obspath, time, var='updraft_helicity', 
        thresh=25., probthresh=0.5, rboxpath=None):
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
                 Variable P_HYD[0,:,:,:]:
                 [0]   Refl > 40 dBZ probs 
                 [1]   UH > 25 m2/s2 probs 
                 [2]   UH > 100 m2/s2 probs
                 [3]   Wind Speed > 40 mph probs
    
    '''
    # Open probabilistic forecast and observational datasets
    probdat = Dataset(probpath)
    obsdat = Dataset(obspath)
    initstr = obsdat.START_DATE
    runinit = datetime(year=int(initstr[:4]), month=int(initstr[5:7]), 
                       day=int(initstr[8:10]))
    fhrs = obsdat.variables['fhr'][:]
    #print(obtimes, time)
    obind = np.where(fhrs-date2num(time, 'hours since ' + initstr) == 0)[0][0]
    print(obind)
    # Choose correct indices based on variable and threshold
    probinds = {'reflectivity' : {40 : 0}, 
                'updraft_helicity' : {25 : 1, 100 : 2},
                'wind_speed' : {40 : 3}}
    # If UH, pull practically perfect
    if var == 'updraft_helicity':
        obs = obsdat.variables['practically_perfect'][:]
        #obs_gt_than = 
    else:
        raise ValueError('Support for {} not yet built in.'.format(var))
    # Pull and splice probability variable
    probvar = probdat.variables['P_HYD'][0]
    d = probinds[var]
    probs = probvar[d[int(thresh)]]
    # First calculate FBS (Fractions Briar Score) on whole grid
    wrf.disable_xarray()
    lats = wrf.getvar(probdat, 'lat')
    lons = wrf.getvar(probdat, 'lon')
    npts = len(lats[:,0]) * len(lons[0,:])f
    #prob_gt_than = (probs >= probthresh)
    fbs = 0.
    fbs_worst = 0.
    print(np.max(probs), np.max(obs[obind]))
    print(np.shape(probs), np.shape(obs[obind]))
    for i in range(len(lons[0,:])):
        for j in range(len(lats[:,0])):
            fbs += (probs[j,i] - obs[obind,j,i])**2
            fbs_worst += probs[j,i]**2 + obs[obind,j,i]**2
            #print(probs[j,i], obs[obind,j,i])
    fbs, fbs_worst = fbs/npts, fbs_worst/npts
    # Use FBS and FBS worst to calculate FSS for whole grid
    fss = 1 - (fbs/fbs_worst)
    print(fss)
    print(npts)
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
        print(rbox)
        #print(lats[mask], lons[mask])
        #print(np.shape(masked_lons))
        print(len(probs[mask]))
        npts = len(probs[mask])
        fbs = 0.
        fbs_worst = 0.
        for i in range(len(masked_probs)):
            fbs += (masked_probs[i] - masked_obs[i])**2
            fbs_worst += masked_probs[i]**2 + masked_obs[i]**2
        fbs, fbs_worst = fbs/npts, fbs_worst/npts
        # Use FBS and FBS worst to calculate FSS for whole grid
        fss = 1 - (fbs/fbs_worst)
    print(fss)    
    #print(rbox)
    return
    
        