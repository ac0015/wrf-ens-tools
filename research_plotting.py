#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Realtime Plotting Functions
---------------------------
To be used with realtime Fortran driver and
the subset/sens objects

Created on Tue Mar 27 10:03:04 2018

@author: aucolema
"""

import matplotlib
matplotlib.use('agg')

import os
import csv
import sys
import numpy as np
from netCDF4 import Dataset
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
#from matplotlib.mlab import griddata
#from interp_analysis import bilinear_interp
import nclcmaps
import matplotlib.patches as patches
from datetime import timedelta, datetime
from copy import copy
import coordinateSystems as cs
import pyproj
import scipy as sp
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage
from subprocess import call
import cmocean
from scipy.ndimage.filters import gaussian_filter

########################################################
# Slice information
# ----------------
# smat/mmat (sensitivity vars and means)
# Smat in P_HYD; Mmat in Q_ICE
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
#
# Rfuncs
#  1 - Avg Sim Refl
#  2 - Max Sim Refl
#  3 - Avg UH
#  4 - Max UH
#  5 - Accum PCP
#  6 - Max Wind Spd
##########################################################

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
    runinitdate: datetime obj for 
        SPC storm reports
    sixhr: boolean specifying whether
        to calculate practically perfect
        probs over six hr time window or 
        use one hr time window
    rtime: time (in num fcst hrs
        from runinit) to obtain six hr
        practically perfect
    Outputs
    -------
    returns tuple containing pract perf probs 
    and lons/lats (respectively) of pperf grid
    '''
    #Get initialization date
    runinitdatef = runinitdate.strftime('%y%m%d')
    rdate = runinitdate + timedelta(hours=rtime)
    if (rtime > 24) & (runinitdate.hour == 12):
        runinitdatef = rdate.strftime('%y%m%d')
    elif (rtime > 35) & (runinitdate.hour == 0):
        runinitdatef = rdate.strftime('%y%m%d')
    #print(runinitdatef)
        
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
    mod = 24
    
    #If there aren't any reports, practically perfect is zero across grid
    if length == 0:
    	pperf = np.zeros_like(wrflon)    
    #Otherwise, let's grid the reports
    else:
        # If six hour, mask reports by valid times in window
        if sixhr:
            hour = rdate.hour
            hours = [(hour - i)%mod for i in range(1,7)]
            mask = [(hr in hours) for hr in time]
        else:
            hour = rdate.hour
            mask = [(hr == (hour-1)%mod) for hr in time]
            print('Reports valid {} to {}'.format((hour-1)%mod,hour))
        try:
            #Set up empty grid onto correct projection
            grid = np.zeros_like(lon)
            NDFD = pyproj.Proj("+proj=lcc +lat_1=25 +lat_2=25 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
        
            # Convert lat/lon grid into projection space
            X, Y = NDFD(lon, lat) 
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
            tmppperf = ndimage.gaussian_filter(grid,sigma=sigma, order=0)
            
            # Interpolate to WRF grid
            pperf = bilinear_interp(X, Y, WRFX, WRFY, tmppperf)
            print('Min/Max Prac Perf: ', np.min(pperf), '/', np.max(pperf))
        except:
            pperf = np.zeros_like(wrflon)
    
    #Remove report CSV file
    os.remove(rptfile)    
    
    return pperf, wrflon, wrflat, np.array(lons)[mask], np.array(lats)[mask]
 

################################################################################################
# Begin plotting modules...
#
# Prob dataset levels (realtime)
#    Lev 1  SLP mean (to add to probability plots for a general synopsis of surface pattern)
#    Lev 2  10-m wind U component mean (to add to probability plots for a general synopsis of surface pattern)
#    Lev 3  10-m wind V component mean (to add to probability plots for a general synopsis of surface pattern)
#    Lev 4  neighborhood prob dBZ exceeds 40 FULL ENSEMBLE
#    Lev 5 - neighborhood prob UH exceeds 25 FULL ENSEMBLE
#    Lev 6 - neighborhood prob UH exceeds 100 FULL ENSEMBLE (we probably won’t use this one)
#    Lev 7 - neighborhood prob sfc wind speed exceeds 40mph FULL ENSEMBLE (we probably won’t use this one)
#    Lev 8  neighborhood prob dBZ exceeds 40 SUBSET based on UH max Response
#    Lev 9 - neighborhood prob UH exceeds 25 SUBSET based on UH max Response
#    Lev 10 - neighborhood prob UH exceeds 100 SUBSET based on UH max Response (we probably won’t use this one)
#    Lev 11 - neighborhood prob sfc wind speed exceeds 40mph SUBSET based on UH max Response (we probably won’t use this one)
#    Lev 12  neighborhood prob dBZ exceeds 40 SUBSET based on dBZ coverage
#    Lev 13 - neighborhood prob UH exceeds 25 SUBSET based on dBZ coverage
#    Lev 14 - neighborhood prob UH exceeds 100 SUBSET based on dBZ coverage (we probably won’t use this one)
#    Lev 15 - neighborhood prob sfc wind speed exceeds 40mph SUBSET based on dBZ coverage (we probably won’t use this one)
#    Lev 16  neighborhood prob dBZ exceeds 40 SUBSET based on UH coverage
#    Lev 17 - neighborhood prob UH exceeds 25 SUBSET based on UH coverage
#    Lev 18 - neighborhood prob UH exceeds 100 SUBSET based on UH coverage (we probably won’t use this one)
#    Lev 19 - neighborhood prob sfc wind speed exceeds 40mph SUBSET based on UH coverage (we probably won’t use this one)    
################################################################################################

def plotPracPerf(runinitdate, sixhr, rtime, sigma=2, outpath='pperf.png'):
    '''
    Plotting SPC practically perfect
    calculations adapted from SPC code
    
    Inputs
    ------
    runinitdate: datetime obj for 
        SPC storm reports
    sixhr: boolean specifying whether
        to calculate practically perfect
        probs over six hr time window or 
        use one hr time window
    rtime: time (in num fcst hrs
        from runinit) to obtain six hr
        practically perfect
    outpath: optional string for 
        absolute path of plot
    Outputs
    -------
    returns NULL and saves plot to outpath
    '''
    pperf, lons, lats, rlons, rlats = calc_prac_perf(runinitdate, sixhr, rtime, sigma=sigma)
    time = runinitdate + timedelta(hours=rtime)
    fig = plt.figure(figsize=(10, 10))
    # Build map
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    ax.set_extent([-108., -85., 28., 45.])
    state_borders = cfeat.NaturalEarthFeature(category='cultural',
               name='admin_1_states_provinces_lakes', scale='50m', facecolor='None') 
    ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
    ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
    ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
    cflevels = np.arange(1, 101, 1)
    cf = ax.contourf(lons, lats, pperf*100., cflevels, transform=ccrs.PlateCarree(),
                     cmap='viridis')
    ax.scatter(rlons, rlats, transform=ccrs.PlateCarree(),
               c='orange', edgecolor='k', alpha=0.8, zorder=10, label='SPC Reports')
    #plt.clabel(cs)
    plt.legend()
    plt.colorbar(cf, ax=ax, fraction=0.039, pad=0.01, orientation='vertical', label='Probability')
    plt.title('Practically Perfect valid: ' + str(time-timedelta(hours=1)) + ' to ' + str(time))
    plt.savefig(outpath)
    return

def plotProbs(probpath, wrfrefpath, rbox, time, nbrhd, outpath='', subset=False):
    '''
    Plots 1-hr probabilities for a specified time for each response
    function (fullens and subset) and overlays the response 
    function box. Only supports 1-hr prob files, as it needs
    the SLP and 10m Wind means to plot.
    
    Inputs
    ------
    probpath ---- string specifying relative or absolute file
                    path of probability netCDF file.
    wrfrefpath -- string specifying relative or absolute file
                    path for an inner domain WRF reference file.
    rbox -------- tuple of floats with resbonse box bounds in the
                    order (llon, ulon, llat, ulat).
    time -------- response time in number of forecast hours (only
                    used for label purposes).
    nbrhd ------- neighborhood of probabilities used (in km)
    outpath ----- string specifying absolute path to store plots.
    '''
    # Get prob data
    probdat = Dataset(probpath)
    probs = probdat.variables['P_HYD'][0]
    refl40 = probs[0]
    uh25 = probs[1]
    uh40= probs[2]
    uh100 = probs[3]
    problist = [refl40, uh25, uh40, uh100]
    
    # Get lat/lon data
    wrfref = Dataset(wrfrefpath)
    clon, clat = wrfref.CEN_LON, wrfref.CEN_LAT
    tlat1, tlat2 = wrfref.TRUELAT1, wrfref.TRUELAT2
    lons, lats = wrfref.variables['XLONG'][0], wrfref.variables['XLAT'][0]
    
    # Build response boxz
    llon, ulon, llat, ulat = rbox
    width = ulon - llon
    height = ulat - llat 

    # Label Strings
    rstrs = ['Reflectivity > 40 dBZ', r'UH > 25 m$^2$/s$^2$', 
             r'UH > 40 m$^2$/s$^2$', r'UH > 100 m$^2$/s$^2$'] 
    if subset:
        figstrs = ['refl40subset', 'uhmax25subset', 'uhmax40subset',
               'uhmax100subset']
    else:
        figstrs = ['refl40fullens', 'uhmax25fullens', 'uhmax40fullens',
               'uhmax100fullens']

    for i in range(len(figstrs)):
        fig = plt.figure(figsize=(10, 10))
        # Build map
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=clon, 
                                                                       central_latitude=clat, 
                                                                       standard_parallels=(tlat1, tlat2)))
        state_borders = cfeat.NaturalEarthFeature(category='cultural',
               name='admin_1_states_provinces_lakes', scale='50m', facecolor='None') 
        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
        # Add rbox and zoom extent to rbox/nearest surrounding area
        rbox = patches.Rectangle((llon, llat), width, height, transform=ccrs.PlateCarree(), 
                             fill=False, color='green', linewidth=2., zorder=3.)
        ax.add_patch(rbox)
        ax.set_extent([llon-10.0, ulon+10.0, llat-5.0, ulat+5.0])
        # Plot probs        
        cflevels = np.linspace(0., 100., 21)
        prob = ax.contourf(lons, lats, gaussian_filter(problist[i], 1), cflevels, transform=ccrs.PlateCarree(), 
                           cmap=nclcmaps.cmap('precip3_16lev'), alpha=0.7, antialiased=True)
        fig.colorbar(prob, fraction=0.046, pad=0.04, orientation='horizontal', label='Probability (Percent)')
        # Format titles and figure names
        ax.set_title(r'Probability of {} at f{} with Neighborhood of {} km'.format(rstrs[i], 
                     time, nbrhd))
        plt.savefig('{}{}nbr{}prob_{}'.format(outpath, figstrs[i], int(nbrhd), time))
        plt.close()
    return
        
def plotDiff(fensprobpath, subprobpath, wrfrefpath, rbox, time,
             responsedate, stormreports=False, outpath=''):
    '''
    Method that plots the difference in 1-hr probabilities between full
    ensemble and subsets for each response function at a specific
    time. (Set up to parse 1-hr prob files only, 6-hr prob files not
    supported.)
    
    Inputs
    ------
    probpath -- string specifying path of full ensemble
                        and subset probability netCDF output file.
    wrfrefpath ------- string specifying path of WRF reference
                        file that contains inner domain data.
    rbox ------------- tuple containing bounds of response function
                        box in order (llon, ulon, llat, ulat).
    responsedate ----- string containing response function date in form
                        2-digit year, month, day (ex. 160508 for May 8
                        2016). Will be used to overlay SPC reports.
    stormreports ----- boolean specifying whether to overlay SPC storm
                        reports for the day.
    '''    
    # Get prob data
    fensprobdat = Dataset(fensprobpath)
    fensprobs = fensprobdat.variables['P_HYD'][0]
    refl40full = fensprobs[0]
    uh25full = fensprobs[1]
    uh40full = fensprobs[2]
    uh100full = fensprobs[3]
    subprobdat = Dataset(subprobpath)
    subprobs = subprobdat.variables['P_HYD'][0]
    refl40sub = subprobs[0]
    uh25sub = subprobs[1]
    uh40sub = subprobs[2]
    uh100sub = subprobs[3]
    fullensprob = [refl40full, uh25full, uh40full, uh100full]
    subsetprob = [refl40sub, uh25sub, uh40sub, uh100sub]
    
    # Get lat/lon data
    wrfref = Dataset(wrfrefpath)
    clon, clat = wrfref.CEN_LON, wrfref.CEN_LAT
    tlat1, tlat2 = wrfref.TRUELAT1, wrfref.TRUELAT2
    lons, lats = wrfref.variables['XLONG'][0], wrfref.variables['XLAT'][0]
    
    # Build response box
    llon, ulon, llat, ulat = rbox
    width = ulon - llon
    height = ulat - llat 
    
    # Specific to fotran program probcalcSUBSETnew. Change if fortran
    #  code changes.
    rstrs = ['Reflectivity > 40 dBZ', r'UH > 25 m$^2$/s$^2$', 
             r'UH > 40 m$^2$/s$^2$', r'UH > 100 m$^2$/s$^2$'] 
    figstrs = ['refl40sub', 'uhmax25sub', 'uhmax40sub', 'uhmax100sub']
    
    # Pull storm reports if needed.
    if stormreports:
        tor = np.atleast_2d(np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                            responsedate+'_rpts_torn.csv', delimiter=',', 
                            skip_header=1, usecols=(5,6), dtype=str))
        hail = np.atleast_2d(np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                            responsedate+'_rpts_hail.csv', delimiter=',', 
                            skip_header=1, usecols=(5,6), dtype=str))
        tlats = tor[:,0]
        tlons = tor[:,1]
        hlats = hail[:,0]
        hlons = hail[:,1]
    
    for i in range(len(rstrs)):
        fig = plt.figure(figsize=(10, 10))
        # Build map projection
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=clon, 
                                                                       central_latitude=clat, 
                                                                       standard_parallels=(tlat1, tlat2)))
        state_borders = cfeat.NaturalEarthFeature(category='cultural',
               name='admin_1_states_provinces_lakes', scale='50m', facecolor='None') 
        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
        # Add response function box
        rbox = patches.Rectangle((llon, llat), width, height, transform=ccrs.PlateCarree(), 
                             fill=False, color='green', linewidth=2., zorder=3.)
        ax.add_patch(rbox)
        ax.set_extent([llon-10.0, ulon+10.0, llat-5.0, ulat+5.0])

        # Set plot difference levels
        cflevels = np.linspace(-80., 80., 161)
        prob = ax.contourf(lons, lats, gaussian_filter(subsetprob[i] - fullensprob[i], 1), cflevels, transform=ccrs.PlateCarree(), 
                           cmap=nclcmaps.cmap('BlWhRe'),alpha=1, antialiased=True)
        if stormreports:
            ax.scatter(np.array(tlons, dtype=float), np.array(tlats, dtype=float),
                       transform=ccrs.PlateCarree(), c='red', edgecolor='k', label='Tor Report', alpha=0.6)
            ax.scatter(np.array(hlons, dtype=float), np.array(hlats, dtype=float), 
                       transform=ccrs.PlateCarree(), c='green', edgecolor='k', label='Hail Report', alpha=0.6)
            plt.legend()
        fig.colorbar(prob, fraction=0.046, pad=0.04, orientation='horizontal', label='Percent Difference')
        ax.set_title(r'Probability difference (Subset - Full Ensemble) of {} at f{}'.format(rstrs[i], str(time)))
        plt.savefig("{}{}probdiff_f{}".format(outpath,figstrs[i], str(time)))
        plt.close()
    return

def plotSPC(outputdir, rbox, responsedate, wrfrefpath, stormreports=False):
    '''
    Plots hourly SPC storm reports up until the response time 
    plus 2 hours (just in case timing of forecast was off).
    
    Inputs
    ------
    outputdir ----- string specifying directory to place hourly
                        storm report output.
    rbox ---------- tuple containing bounds of response function
                        box in order (llon, ulon, llat, ulat).
    responsedate -- string containing response function date in form
                        2-digit year, month, day (ex. 160508 for May 8
                        2016). Will be used to overlay SPC reports.
    wrfrefpath ---- string specifying file path of WRF file that
                        contains lat/lon info on inner domain
    ''' 
    try:
        # Pull tornado and hail report csvs, only taking necessary columns
        tor = np.atleast_2d(np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                                str(responsedate)+'_rpts_torn.csv',
                                delimiter=',', skip_header=1, usecols=(5,6), dtype=str))
        hail = np.atleast_2d(np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                                str(responsedate)+'_rpts_hail.csv', 
                                delimiter=',', skip_header=1, usecols=(5,6), dtype=str))
        # Splice csv files into locs, and times
        torlats = tor[:,0]
        torlons = tor[:,1]
        hlats = hail[:,0]
        hlons = hail[:,1]
        
        # Get lat/lon data for plot extent
        wrfref = Dataset(wrfrefpath)
        clon, clat = wrfref.CEN_LON, wrfref.CEN_LAT
        tlat1, tlat2 = wrfref.TRUELAT1, wrfref.TRUELAT2
        
        # Build response box
        llon, ulon, llat, ulat = rbox
        width = ulon - llon
        height = ulat - llat 
        
        # Plot
        fig = plt.figure(figsize=(10, 10))
        
        # Build projection/map
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=clon, 
                                                                           central_latitude=clat, 
                                                                           standard_parallels=(tlat1, tlat2)))
        state_borders = cfeat.NaturalEarthFeature(category='cultural',
                   name='admin_1_states_provinces_lakes', scale='50m', facecolor='None') 
        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
        # Add rbox and zoom extent to rbox/nearest surrounding area
        rbox = patches.Rectangle((llon, llat), width, height, transform=ccrs.PlateCarree(), 
                                 fill=False, color='green', linewidth=2., zorder=3.)
        ax.add_patch(rbox)
        ax.set_extent([llon-10.0, ulon+10.0, llat-5.0, ulat+5.0])
        
        # Add reports as scatter points
        ax.scatter(np.array(torlons, dtype=float), np.array(torlats, dtype=float), 
                   transform=ccrs.PlateCarree(), c='red', edgecolor='k', 
                   label='Tor Report', alpha=0.3)
        ax.scatter(np.array(hlons, dtype=float), np.array(hlats, dtype=float), 
                   transform=ccrs.PlateCarree(), c='green', edgecolor='k', 
                   label='Hail Report', alpha=0.3)
        ax.set_title(r'SPC Tor and Hail Reports valid {}'.format(str(responsedate)))
        plt.savefig("{}SPCreport{}".format(outputdir, str(responsedate)))
        plt.close()
    except:
        print("No hail or tornado reports to plot")
    return

def calc1hrPaintball(datem, hour, variable, members, nx, ny, thresh):
    hour='%02d' % int(hour)
    start=int(hour)-1; start='%02d' % start

    numens = len(members)
    memcount = 0
    base='/lustre/research/bancell/aucolema/HWT2016runs/' + datem  +'/'
    
    ntimes=[1]

    ####################################################################
    max_mem_val=np.empty(numens)
    paintball=np.empty((numens,ny,nx),dtype='float') #allocate numens x model grid array
    for i in range(len(members)): ###loop through num of members      
        path = base + 'mem' + str(members[i]) + '/' + 'R' + str(members[i]) + '_' + hour + '.out'
        if os.path.isfile(path):
            var=np.empty((len(ntimes),ny,nx),dtype='float')
    
            memcount += 1
            #print('check ' + str(i))
            
            #### Time Loop, num of wrfouts ####
            for j in ntimes: ###loop through num of files/times
    
                k=ntimes.index(j) #loop count    
                ncfile = Dataset(path, 'r')
                if variable.lower() == 'refl_10cm':
                    var[k,:,:] = ncfile.variables[variable][0,0,:,:] #for 4d refl var
                elif variable.lower() == 'up_heli_max': 
                    var[k,:,:] = ncfile.variables[variable][0,:,:] #for 3d uh var
                elif variable.lower() == 'wpsd10max':
                    var[k,:,:] = ncfile.variables[variable][0,:,:] #for 3d wspd var
                ncfile.close()
      ##################################
    
        ### Back to Each Member ###
            #Now we have our NT x NY x NX array for member i...find max values!
            max=np.max(var,0) #This is max 2d field over three times
            
            max_mem_val[i]=np.max(np.max(max,0))
            #print(max_mem_val[i])    
            
            #give occurence integer of mem number
            max[max > thresh] = int(members[i]) #if >= thresh, set as 1
            max[max != int(members[i])] = 0 # if not 1 (set from thresh), set as 0
            paintball[i,:,:]=max
            
    return paintball

def calc6hrPaintball(datem, hour, variable, members, nx, ny, varstr):
    hour='%02d' % int(hour)
    start=int(hour)-6; start='%02d' % start

    if varstr == 'UH':
         thresh=25 #threshold grid cells will be paintballed on map
    elif varstr == 'DBZ':
         thresh=40
    numens = len(members)
    memcount = 0
    base='/home/bancell/RT_ENS/holdver/' + datem + '/'
    print(base)

    ######## DATETIME PROCESSING ########################################
    h1=int(hour)-5; h2=int(hour)-4 ; h3=int(hour)-3; h4=int(hour)-2; 
    h5=int(hour)-1; h6=int(hour);
    
    time1 = datetime.strptime(datem, '%Y%m%d%H')
    time1 += timedelta(hours=h1)
    time1=time1.strftime('%Y%m%d%H')
    
    time2 = datetime.strptime(datem, '%Y%m%d%H')
    time2 += timedelta(hours=h2)
    time2=time2.strftime('%Y%m%d%H')
    
    time3 = datetime.strptime(datem, '%Y%m%d%H')
    time3 += timedelta(hours=h3)
    time3=time3.strftime('%Y%m%d%H')
    
    time4 = datetime.strptime(datem, '%Y%m%d%H')
    time4 += timedelta(hours=h4)
    time4=time4.strftime('%Y%m%d%H')
    
    time5 = datetime.strptime(datem, '%Y%m%d%H')
    time5 += timedelta(hours=h5)
    time5=time5.strftime('%Y%m%d%H')
    
    time6 = datetime.strptime(datem, '%Y%m%d%H')
    time6 += timedelta(hours=h6)
    time6=time6.strftime('%Y%m%d%H')
    
    yyyy1=time1[0:4]; mm1=time1[4:6]; dd1=time1[6:8]; hh1=time1[8:10]
    yyyy2=time2[0:4]; mm2=time2[4:6]; dd2=time2[6:8]; hh2=time2[8:10]
    yyyy3=time3[0:4]; mm3=time3[4:6]; dd3=time3[6:8]; hh3=time3[8:10]
    yyyy4=time4[0:4]; mm4=time4[4:6]; dd4=time4[6:8]; hh4=time4[8:10]
    yyyy5=time5[0:4]; mm5=time5[4:6]; dd5=time5[6:8]; hh5=time5[8:10]
    yyyy6=time6[0:4]; mm6=time6[4:6]; dd6=time6[6:8]; hh6=time6[8:10]
    
    wrf1='wrfout_d02_'+ str(yyyy1) + '-' + str(mm1) + '-' + str(dd1) + '_' + str(hh1) + ':00:00'
    wrf2='wrfout_d02_'+ str(yyyy2) + '-' + str(mm2) + '-' + str(dd2) + '_' + str(hh2) + ':00:00'
    wrf3='wrfout_d02_'+ str(yyyy3) + '-' + str(mm3) + '-' + str(dd3) + '_' + str(hh3) + ':00:00'
    wrf4='wrfout_d02_'+ str(yyyy4) + '-' + str(mm4) + '-' + str(dd4) + '_' + str(hh4) + ':00:00'
    wrf5='wrfout_d02_'+ str(yyyy5) + '-' + str(mm5) + '-' + str(dd5) + '_' + str(hh5) + ':00:00'
    wrf6='wrfout_d02_'+ str(yyyy6) + '-' + str(mm6) + '-' + str(dd6) + '_' + str(hh6) + ':00:00'
    
    a=[wrf1,wrf2,wrf3,wrf4,wrf5,wrf6]
    ####################################################################
    max_mem_val=np.empty(numens)
    paintball=np.empty((numens,ny,nx),dtype='float') #allocate numens x model grid array
    for i in range(len(members)): ###loop through num of members      
        if os.path.isfile(base + 'mem' + str(members[i]) + '/' + a[-1]):
            var=np.empty((len(a),ny,nx),dtype='float')
    
            memcount += 1
            
            #### Time Loop, num of wrfouts ####
            for j in a: ###loop through num of files/times
    
                k=a.index(j) #loop count
                path=base+'mem'+str(members[i])+"/"+j
                #print(path)
    
                ncfile = Dataset(path, 'r')
                if variable.lower() == 'refl_10cm':
                    var[k,:,:] = ncfile.variables[variable][0,0,:,:] #for 4d refl var
                elif variable.lower() == 'up_heli_max': 
                    var[k,:,:] = ncfile.variables[variable][0,:,:] #for 3d uh var
                ncfile.close()
      ##################################
    
        ### Back to Each Member ###
            #Now we have our NT x NY x NX array for member i...find max values!
            max=np.max(var,0) #This is max 2d field over three times
            
            max_mem_val[i]=np.max(np.max(max,0))
            #print("Member " + str(i) + " Max:", max_mem_val[i])    
            
            #give occurence integer of mem number
            max[max > thresh] = int(members[i]) #if >= thresh, set as 1
            max[max != int(members[i])] = 0 # if not 1 (set from thresh), set as 0
            paintball[i,:,:]=max
            
    return paintball

def plotHrlySPC(outputdir, runinit, rbox, numtimes, wrfrefpath):
    '''
    Plots hourly SPC storm reports up until the response time 
    plus 2 hours (just in case timing of forecast was off). 
    Note - only plots tornado and hail right now. More
    for verification with UH.
    
    Inputs
    ------
    outputdir ----- string specifying directory to place hourly
                        storm report output.
    runinit ------- datetime object for the model run initialization
                        time (start of SPC hourly)
    numtimes ------ integer for number of forecast hours to plot
                        (max is 24, only for plotting one day's
                        worth of SPC storm reports)
    wrfrefpath ---- string specifying file path of WRF file that
                        contains lat/lon info on inner domain
    '''
    # Slice time object
    yr, mo, day, hr, mn, sec, wday, yday, isdst = runinit.timetuple()
    if len(str(mo)) < 2: mo = '0' + str(mo)
    if len(str(day)) < 2: day = '0' + str(day)
    
    try:
        # Pull tornado and hail report csvs, only taking necessary columns
        tor = np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                                str(yr)[-2:]+str(mo)+str(day)+'_rpts_torn.csv',
                                delimiter=',', skip_header=1, usecols=(0,5,6), dtype=str)
        hail = np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                                str(yr)[-2:]+str(mo)+str(day)+'_rpts_hail.csv', 
                                delimiter=',', skip_header=1, usecols=(0,5,6), dtype=str)
        
        # Splice csv files into locs, and times
        torlats = tor[:,1]
        torlons = tor[:,2]
        tortimes = tor[:,0]
        hlats = hail[:,1]
        hlons = hail[:,2]
        hailtimes = hail[:,0]
        
        # Splice times by hour
        torhr = []
        hailhr = []
        for tortime in tortimes: torhr.append(tortime[:2])
        for hailtime in hailtimes: hailhr.append(hailtime[:2])
        torhr = np.array(torhr, dtype=int)
        hailhr = np.array(hailhr, dtype=int)
        
        # Create list of times to plot
        maxtimes = 24
        if numtimes > maxtimes: 
            numtimes = 24 # Limit number of hours to full day's worth of reports
        times = np.arange(hr, hr+numtimes) % 24
        
        # Get lat/lon data for plot extent
        wrfref = Dataset(wrfrefpath)
        clon, clat = wrfref.CEN_LON, wrfref.CEN_LAT
        tlat1, tlat2 = wrfref.TRUELAT1, wrfref.TRUELAT2
        
        # Build response box
        llon, ulon, llat, ulat = rbox
        width = ulon - llon
        height = ulat - llat 
    
        # Plot
        for i in range(len(times)):
            # Format like SPC (reports run from 12Z day 1 to 1159Z day 2)
            if times[i] < 12: 
                date = runinit + timedelta(days=1, hours=int(times[i]))
            else:
                date = runinit + timedelta(hours=int(times[i]))
            yr, mo, day, hr, mn, sec, wday, yday, isdst = date.timetuple()
            # Match times with masks
            tmask = (torhr[:] == times[i])
            hmask = (hailhr[:] == times[i])
            # Plot
            fig = plt.figure(figsize=(10, 10))
            # Build projection/map
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(central_longitude=clon, 
                                                                           central_latitude=clat, 
                                                                           standard_parallels=(tlat1, tlat2)))
            state_borders = cfeat.NaturalEarthFeature(category='cultural',
                   name='admin_1_states_provinces_lakes', scale='50m', facecolor='None') 
            ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
            ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
            ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
            # Add rbox and zoom extent to rbox/nearest surrounding area
            rbox = patches.Rectangle((llon, llat), width, height, transform=ccrs.PlateCarree(), 
                                 fill=False, color='green', linewidth=2., zorder=3.)
            ax.add_patch(rbox)
            ax.set_extent([llon-10.0, ulon+10.0, llat-5.0, ulat+5.0])
            # Add tor/hail pts
            if (len(torlons[tmask]) > 0):
                ax.scatter(np.array(torlons[tmask], dtype=float), np.array(torlats[tmask], dtype=float), 
                           transform=ccrs.PlateCarree(), c='red', edgecolor='k', label='Tor Report', alpha=0.7,
                           zorder=12)
            if (len(hlons[hmask]) > 0):
                ax.scatter(np.array(hlons[hmask], dtype=float), np.array(hlats[hmask], dtype=float), 
                           transform=ccrs.PlateCarree(), c='green', edgecolor='k', label='Hail Report', alpha=0.7)
            plt.legend()
            # Name and save
            if len(str(mo)) < 2: mo = '0' + str(mo)
            if len(str(day)) < 2: day = '0' + str(day)
            if len(str(times[i])) < 2: timestr ='0' + str(times[i])
            else: timestr = str(times[i])
            ax.set_title(r'SPC Storm Reports valid {} UTC to {} UTC'.format(str(date+timedelta(hours=-1)), str(date)))
            plt.savefig("{}SPCreport{}{}{}_{}Z".format(outputdir, str(yr), str(mo), 
                        str(day), timestr))
            plt.close()
    except:
        print("No reports to plot...")
    return
 
def plotSixPanels(dirdate, stormreports, submems, sixhour=True, time=None, 
                  subsettype='uhmax', nbrhd=30):
    '''
    Plots three six panel plots for three response functions
    at a specified time and area around a response box in the following 
    order:
        (1,1) Full Ensemble Probabilities
        (1,2) Subset Probabilities
        (2,1) Probability Differences
        (2,2) SPC Reports/Practically Perfect Probs
        (3,1) Full Ensemble Paintball
        (3,2) Subset Paintball
    
    Inputs
    ------
    dirdate ------ model run initialization string (YYYYMMDD)
    rbox --------- tuple of floats with resbonse box bounds in the
                    order (llon, ulon, llat, ulat).
    stormerports - boolean specifying whether to plot storm reports
                    from SPC or not
    submems --- list with members to use for subset
    sixhour ------ boolean specifying whether using the six hour probs
                    and paintball (False is one hour). Defaults to True
    time --------- optional integer specifying time to plot (in forecast
                    hours). If left None, uses rtime from esens.in
    subsettype --- optional string specifying response function used to 
                    subset. Only affects figure namefor organizational purposes.
    '''    
    # Base dir
    yr, mo, day, hr = str(dirdate)[:4], str(dirdate)[4:6], str(dirdate)[6:8], str(dirdate)[8:10]
    runinit = datetime(year=int(yr), month=int(mo), day=int(day), hour=int(hr))
    base='/lustre/research/bancell/aucolema/HWT2016runs/' 
    
    # Get subset data
    subsetdat = np.genfromtxt(base + dirdate + '/esens.in', dtype=str)
    # If time param not overidden, use time from subset data
    if time == None:
        time = int(subsetdat[2])
    rbox = np.array(subsetdat[4:8], dtype=float)
    print("RTime: ", time, "RBox: ", rbox)
    
    # TO-DO: Replace generic probppath with 1-hr and 6-hr paths/probs
    if sixhour:
        fullensprobpath = base + dirdate + '/probs/FULLENSwrfout_nbr{}_f{}.prob'.format(str(int(nbrhd)), str(time))
    else:
        fullensprobpath = base + dirdate + '/probs/FULLENSwrfout_nbr{}_f{}.prob'.format(str(int(nbrhd)),str(time))
        subsetprobpath = base + dirdate + '/probs/SUBSETwrfout_nbr{}_f{}.prob'.format(str(int(nbrhd)),str(time))
    wrfrefpath = base + dirdate + '/wrfoutREFd2'
    
    # Get prob data
    # Six hour doesn't have means, so must go to one hour data for it
    # TO-DO: Rewrite for research sixhour - started but unfinished
    if sixhour:
        meanpath = base + dirdate + '/Rmean.out'
        meandat = Dataset(meanpath)
        onehrprobs = meandat.variables['P_HYD'][0]
        slpmean = onehrprobs[0]/100.
        u10mean = onehrprobs[1]
        v10mean = onehrprobs[2]
        refl40fullens = probs[0]
        uh25fullens = probs[1]
        refl40subuhmax = probs[7]
        uh25subuhmax = probs[8]
        refl40subdbzcov = probs[11]
        uh25subdbzcov = probs[12]
        refl40subuhcov = probs[15]
        uh25subuhcov = probs[16]
        timeframe = '6 hr'
    # Else pull all data from one hour probs
    else:
        meanpath = base + dirdate + '/Rmean.out'
        meandat = Dataset(meanpath)
        sfcpmean_unmasked = meandat.variables['PSFC'][0,:,:]/100.
        u10mean_unmasked = meandat.variables['U10'][0,:,:]
        v10mean_unmasked = meandat.variables['V10'][0,:,:]
        sfcpmean = np.ma.masked_array(sfcpmean_unmasked, mask=(sfcpmean_unmasked>=9e9))
        u10mean = np.ma.masked_array(u10mean_unmasked, mask=(u10mean_unmasked>=9e9))
        v10mean = np.ma.masked_array(v10mean_unmasked, mask=(v10mean_unmasked>=9e9))
        fullensprobdat = Dataset(fullensprobpath)
        fensprobvar = fullensprobdat.variables['P_HYD'][0,:,:,:]
        subsetprobdat = Dataset(subsetprobpath)
        subsetprobvar = subsetprobdat.variables['P_HYD'][0,:,:,:]
        refl40fullens = fensprobvar[0]
        uh25fullens = fensprobvar[1]
        uh40fullens = fensprobvar[2]
        uh100fullens = fensprobvar[3]
        wspd40fullens = fensprobvar[4]
        refl40sub = subsetprobvar[0]
        uh25sub = subsetprobvar[1]
        uh40sub = subsetprobvar[2]
        uh100sub = subsetprobvar[3]
        wspd40sub = subsetprobvar[4]
        timeframe = '1 hr'

    fullensprobs = [uh25fullens, refl40fullens, uh40fullens, uh100fullens, wspd40fullens]
    subsetprobs = [uh25sub, refl40sub, uh40sub, uh100sub, wspd40sub]

    # Create time to plot
    rdate = runinit + timedelta(hours=time)
    hour = rdate.hour
    
    if stormreports:
        # Pull tornado and hail report csvs, only taking necessary columns
        tor = np.atleast_2d(np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                                str(yr)[-2:]+str(mo)+str(day)+'_rpts_torn.csv',
                                delimiter=',', skip_header=1, usecols=(0,5,6), dtype=str))
        hail = np.atleast_2d(np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                                str(yr)[-2:]+str(mo)+str(day)+'_rpts_hail.csv', 
                                delimiter=',', skip_header=1, usecols=(0,1,5,6), dtype=str))
        wind = np.atleast_2d(np.genfromtxt('http://www.spc.noaa.gov/climo/reports/'+
                                str(yr)[-2:]+str(mo)+str(day)+'_rpts_wind.csv', 
                                delimiter=',', skip_header=1, usecols=(0,1,5,6), dtype=str))
        torempty = (not bool(tor.size)); hailempty = (not bool(hail.size)); windempty = (not bool(wind.size))
            
        # Format like SPC (reports run from 12Z day 1 to 1159Z day 2)
        if hour < 12: 
            date = rdate + timedelta(days=-1)
        else:
            date = rdate
        spcyr, spcmo, spcday, spchr, mn, sec, wday, yday, isdst = date.timetuple()
            
        if not torempty:
            # Splice csv files into locs, and times
            torlats = tor[:,1]
            torlons = tor[:,2]
            tortimes = tor[:,0] # Splice times by hour
            torhr = []
            for tortime in tortimes: torhr.append(tortime[:2])
            torhr = np.array(torhr, dtype=int)
            if sixhour:
                # Match six-hr times with masks
                hours = [(hour - i)%24 for i in range(1,7)]
                print("Timeframe used for sixhrlys:", hours)
                tmask = [(hr in hours) for hr in torhr]
            else:
                # Match times with masks
                tmask = (torhr[:] == (hour-1)%24)
        if not hailempty:
            hailtimes = hail[:,0]
            hsizes = np.array(hail[:,1], dtype=float)
            hlats = hail[:,2]
            hlons = hail[:,3]
            hailhr = []
            for hailtime in hailtimes: hailhr.append(hailtime[:2])
            hailhr = np.array(hailhr, dtype=int)
            if sixhour:
                # Match six-hr times with masks
                hours = [(hour - i)%24 for i in range(1,7)]
                hmask = [(hr in hours) for hr in hailhr]
            else:
                hmask = (hailhr[:] == (hour-1)%24)
        if not windempty:
            windtimes = wind[:,0]
            wspd = wind[:,1]
            wlats = wind[:,2]
            wlons = wind[:,3]
            windhr = []
            for windtime in windtimes: windhr.append(windtime[:2])
            windhr = np.array(windhr, dtype=int)
            if sixhour:
                # Match six-hr times with masks
                hours = [(hour - i)%24 for i in range(1,7)]
                wmask = [(hr in hours) for hr in windhr]
            else:
                wmask = (windhr[:] == (hour-1)%24)

        # Create hail size/wind masks to differentiate large hail and high wind
        if not hailempty:
            # Create hail size mask (to plot large hail differently)
            hsizemask = (hsizes > 200.)
        if not windempty:
            # Create wind speed mask (to plot high winds differently)
            wspdmask = np.empty(np.shape(wspd), dtype=bool)
            for i in range(len(wspd)):
                if (wspd[i] != 'UNK'):
                    # Wind is in miles per hour, not knots
                    # >65 kts (>74.8 mph) = high wind
                    if float(wspd[i]) >= 74.8:
                        wspdmask[i] = True
                    else:
                        wspdmask[i] = False
                else:
                    wspdmask[i] = False
            
    # Get lat/lon data
    wrfref = Dataset(wrfrefpath)
    clon, clat = wrfref.CEN_LON, wrfref.CEN_LAT
    tlat1, tlat2 = wrfref.TRUELAT1, wrfref.TRUELAT2
    lons, lats = wrfref.variables['XLONG'][0], wrfref.variables['XLAT'][0]
    wrfref.close()
    
    # Build response box
    llon, ulon, llat, ulat = rbox
    width = ulon - llon
    height = ulat - llat
    
    # Get subset members
    memslist = [submems, submems, submems, submems]

    # Label Strings
    wrfvar = ['UP_HELI_MAX', 'REFL_10CM', 'UP_HELI_MAX', 'UP_HELI_MAX', 'WSPD10MAX']
    thresh = [25, 40, 40, 100, 40]
    titlevars = [r'Updraft Helicity Maximum', 'Reflectivity Average', 
                 r'Updraft Helicity Maximum', r'Updraft Helicity Maximum', 'Wind Speed Maximum']
    rstrs = [r'UH > 25 m$^2$/s$^2$', 'Reflectivity > 40 dBZ', r'UH > 40 m$^2$/s$^2$'
             r'UH > 100 m$^2$/s$^2$', r'Wind Speed > 40 miles/hour'] 
    figstrs = ['uhmax25_sub{}_{}_f{}'.format(subsettype, dirdate, time), 
               'refl40_sub{}_{}_f{}'.format(subsettype, dirdate, time),
               'uhmax40_sub{}_{}_f{}'.format(subsettype, dirdate, time),
               'uhmax100_sub{}_{}_f{}'.format(subsettype, dirdate, time), 
               'wspd40_sub{}_{}_f{}'.format(subsettype, dirdate, time)]
    paintballstrs = ['UH', 'DBZ', 'UH', 'UH', 'Wind Speed']
  
    if sixhour:
        figstrs = ['sixhr' + x for x in figstrs]

    # Pull full ensemble mean and subset ensemble means
    rvals = Dataset(base + dirdate + '/Rvals.nc')
    submems = np.array(submems)
    fullensmeanuh = np.mean(np.array(rvals.variables['UH_MAX'][:]))
    submeanuh = np.mean(np.array(rvals.variables['UH_MAX'][submems-1]))
    fullensmeandbz = np.mean(np.array(rvals.variables['DBZ_MAX'][:]))
    submeandbz = np.mean(np.array(rvals.variables['DBZ_MAX'][submems-1]))
    fullensmeanwind = np.mean(np.array(rvals.variables['WSPD_AVG'][:]))
    submeanwind = np.mean(np.array(rvals.variables['WSPD_AVG'][submems-1]))
    fullmeans = [fullensmeanuh, fullensmeandbz, fullensmeanuh, fullensmeanuh, fullensmeanwind]
    submeans = [submeanuh, submeandbz, submeanuh, submeanuh, submeanwind]
    meandescr = ['UH Max', 'DBZ Max', 'UH Max', 'UH Max', 'Wind Speed Avg']

    # Calculate SPC Practically Perfect Probs for plotting
    pperf, plons, plats = calc_prac_perf(runinit, sixhour, time)

    # Begin plotting madness
    for i in range(len(figstrs)):
        print("Plotting {}".format(figstrs[i]))
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row',
             figsize=(10,12), subplot_kw={'projection': ccrs.LambertConformal(central_longitude=clon, central_latitude=clat, 
                                          standard_parallels=(tlat1, tlat2))})
        axes = [ax1, ax2, ax3, ax4, ax5, ax6]
        mems = np.array(memslist[i][:], dtype=int)
        inds = mems.argsort()
        print('Ordered subset members:', mems[inds])
        cflevs= [0., 2., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.]

        # Add background data/set extents for each plot
        for ax in axes:
            state_borders = cfeat.NaturalEarthFeature(category='cultural',
               name='admin_1_states_provinces_lakes', scale='50m', facecolor='None')
            ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
            ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
            ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
            # Add rbox and zoom extent to rbox/nearest surrounding area
            rbox = patches.Rectangle((llon, llat), width, height, transform=ccrs.PlateCarree(), 
                                 fill=False, color='green', linewidth=2., zorder=3.)
            ax.add_patch(rbox)
            ax.set_extent([llon-6., ulon+6., llat-3., ulat+3.])
            if (ax != ax3) and (ax != ax4) and (ax != ax5) and (ax != ax6):
                # Set clevels for means
                #sfcplev = np.arange(970, 1051, 4)
                # Plot means
                #print(np.min(sfcpmean), np.max(sfcpmean))
                #sfcp = ax.contour(lons, lats, sfcpmean, sfcplev, transform=ccrs.PlateCarree(), colors='k')
                #plt.clabel(sfcp, fmt='%i')
                ax.barbs(lons[::25,::25], lats[::25,::25], 
                          u10mean[::25,::25], v10mean[::25,::25], length=5, 
                          linewidth=0.5, transform=ccrs.PlateCarree(), zorder=4)
        if stormreports & (time>11):
            if not hailempty:
                # Need to isolate both time and hailsize, so combine masks
                reghailmask = hmask & (hsizemask == False)
                lghailmask =  hmask & hsizemask
            if not windempty:
                regwspdmask = wmask & (wspdmask == False)
                highwspdmask = wmask & wspdmask
            # Add reports as scatter points, filtered by time
            if not torempty:
                if (len(torlons[tmask]) > 0):
                    ax4.scatter(np.array(torlons[tmask], dtype=float), np.array(torlats[tmask], dtype=float), 
                               transform=ccrs.PlateCarree(), c='red', edgecolor='k', 
                               label='Tor Report', alpha=0.8, zorder=6)
            if not hailempty:
                if (len(hlons[reghailmask]) > 0):
                    ax4.scatter(np.array(hlons[reghailmask], dtype=float), np.array(hlats[reghailmask], dtype=float), 
                           transform=ccrs.PlateCarree(), c='green', edgecolor='k', 
                           label='Hail Report', alpha=0.8, zorder=4)
                if (len(hlons[lghailmask]) >0):
                    ax4.scatter(np.array(hlons[lghailmask], dtype=float), np.array(hlats[lghailmask], dtype=float), 
                          transform=ccrs.PlateCarree(), c='k', marker='^',
                          edgecolor='k', label='Large Hail Report (>2")', alpha=0.8, zorder=5)
            if not windempty:
                if (len(wlons[regwspdmask]) > 0):
                    ax4.scatter(np.array(wlons[regwspdmask], dtype=float), np.array(wlats[regwspdmask], dtype=float), 
                           transform=ccrs.PlateCarree(), c='blue', edgecolor='k', 
                           label='Wind Report', alpha=0.8, zorder=2)
                if (len(wlons[highwspdmask]) > 0):
                    ax4.scatter(np.array(wlons[highwspdmask], dtype=float), np.array(wlats[highwspdmask], dtype=float), 
                          transform=ccrs.PlateCarree(), c='k', marker="s", edgecolor='k', 
                          label='High Wind Report (>65 kts)', alpha=0.8, zorder=3)
            # Plot practically perfect if values above zero
            if np.max(pperf) > 0.:
                print(np.max(pperf*100))
                cfpperf = ax4.contourf(plons, plats, pperf*100, cflevs, cmap=nclcmaps.cmap('precip3_16lev'), 
                                       transform=ccrs.PlateCarree(), alpha=0.8, zorder=1)
                pperfcbar = fig.colorbar(cfpperf, fraction=0.046, pad=0.04, ax=ax4, orientation='vertical', 
                     label='Probability (Percent)')
                #pperfcbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in cflevs])
            # Add legend for SPC reports
            leg = ax4.legend(loc='upper right', prop={'size': 5})
            leg.set_zorder(20)
            if sixhour:
                ax4.set_title('SPC Reports and Practically Perfect Probs \n Valid {} to {}'.format(str(rdate-timedelta(hours=6)),
                          rdate),fontsize=9)
            else:
                ax4.set_title('SPC Reports and Practically Perfect Probs \n Valid {} to {}'.format(str(rdate-timedelta(hours=1)),
                          rdate),fontsize=9)

        else:
            ax4.set_title('SPC Reports (currently unavailable)')
        print("Max Full Ens Probs: ", np.max(fullensprobs[i]), "Max Subset Probs: ", np.max(subsetprobs[i]))
        # Plot probs        
        fullprob = ax1.contourf(lons, lats, fullensprobs[i], cflevs, transform=ccrs.PlateCarree(), 
                           cmap=nclcmaps.cmap('precip3_16lev'), zorder=1, antialiased=True)
        subprob = ax2.contourf(lons, lats, subsetprobs[i], cflevs, transform=ccrs.PlateCarree(), 
                           cmap=nclcmaps.cmap('precip3_16lev'), zorder=1, antialiased=True)
        deltalevs = np.linspace(-80., 80., 17)
        deltacmap = copy(nclcmaps.cmap('ViBlGrWhYeOrRe'))
        deltaprob = ax3.contourf(lons, lats, (subsetprobs[i] - fullensprobs[i]), 
                                 deltalevs, transform=ccrs.PlateCarree(), 
                                 cmap=deltacmap, 
                                 antialiased=True)
        # Plot paintball
        try:
            fullensrange = np.arange(1,43,1)
            if sixhour:
                fullpaintball = calc6hrPaintball(dirdate, time, wrfvar[i], fullensrange, 
                                         len(lons[0,:]), len(lats[:,0]), paintballstrs[i])
                subpaintball = calc6hrPaintball(dirdate, time, wrfvar[i],mems[inds], 
                                         len(lons[0,:]), len(lats[:,0]), paintballstrs[i])
            else:
                fullpaintball = calc1hrPaintball(dirdate, time, wrfvar[i], fullensrange, 
                                         len(lons[0,:]), len(lats[:,0]), thresh[i])
                subpaintball = calc1hrPaintball(dirdate, time, wrfvar[i],mems[inds], 
                                         len(lons[0,:]), len(lats[:,0]), thresh[i])
            flevs = np.arange(0.99999, len(fullensrange)+0.01,1)
            sublevs = mems[inds]
            if len(sublevs) < 2:
                raise
            for j in range(len(fullensrange)):
                if np.max(fullpaintball[j]) >= flevs[0]: 
                    fpaint = ax5.contourf(lons, lats, fullpaintball[j],
                                        levels=flevs, cmap=plt.cm.jet,
                                        transform=ccrs.PlateCarree(), 
                                        antialiased=True, alpha=0.7)
                    #fpaint2 = ax5.contour(lons, lats, fullpaintball[j],
                    #                    levels=flevs, transform=ccrs.PlateCarree())
                    #inds = np.where(fullpaintball[j] > 0.)
                    #fpaint = ax5.scatter(lons[inds], lats[inds], alpha=0.3, cmap=plt.cm.jet,
                    #                    transform=ccrs.PlateCarree(), s=10.)
            for j in range(len(submems)):
                if np.max(subpaintball[j]) > sublevs[0]:
                    print(np.max(subpaintball[j]), j) 
                    spaint = ax6.contourf(lons, lats, subpaintball[j],
                                        levels=sublevs, cmap=plt.cm.jet,
                                        transform=ccrs.PlateCarree(), 
                                        antialiased=True, alpha=0.7)
            if np.max(fullpaintball) > flevs[0]:
                fpaintballcbar = fig.colorbar(fpaint, fraction=0.046, pad=0.04, orientation='vertical',
                                     ticks=np.arange(0.99999, len(fullensrange)+0.01, 6), 
                                     ax=ax5,label='Full Ens Member Number')
                fpaintballcbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(1,len(fullensrange)+0.05,6)])
            if np.max(subpaintball) > sublevs[0]:
                spaintballcbar = fig.colorbar(spaint, fraction=0.046, pad=0.04, orientation='vertical',
                                          ax=ax6,label='Subset Member Number')
        except:
            print('No paintball vals to plot. Skipping')
        fig.colorbar(fullprob, fraction=0.046, pad=0.04, ax=ax2, orientation='vertical', 
                     label='Probability (Percent)')
        fig.colorbar(deltaprob, fraction=0.046, pad=0.04, ax=ax3, orientation='vertical', 
                     label='Probability (Percent)', extend='both')
        fig.suptitle('Response Function: {} {} at f{} \n Valid for Run Initialized: {} \n Response Endtime: {}'.format(timeframe, titlevars[i], 
                  time, str(runinit), str(rdate)))
        ax1.set_title('Full Ens Prob of {} with {} km Neighborhood \n Mean {}: {:.2f}'.format(rstrs[i],
                      int(nbrhd), meandescr[i], fullmeans[i]))
        ax2.set_title('Subset Prob of {} with {} km Neighborhood \n Mean {}: {:.2f}'.format(rstrs[i], 
                      int(nbrhd), meandescr[i], submeans[i]))
        ax3.set_title(r'Delta Probs (Subset - Full Ensemble)')
        ax5.set_title('Full Ens Paintball of {}'.format(rstrs[i]))
        ax6.set_title('Subset Paintball of {}'.format(rstrs[i]))
        plt.savefig(base + dirdate + '/sixpanel_{}'.format(figstrs[i]))
        plt.close()
    
    return
           
def plot1hrSixPanels(dirdate, stormreports=False, numtimes=48):
    '''
    Plot 1 hrly six panels for three response functions as
    denoted in plotSixPanels()
    Inputs
    -------
    dirdate ------ model run initialization string (YYYYMMDD)
    numtimes ----- number of frames (defaults to traditional 48 fhrs)
    '''
    for i in range(numtimes):
        plotSixPanels(dirdate, stormreports=stormreports, sixhour=False, time=i)

