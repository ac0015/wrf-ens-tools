#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 16:07:43 2018

A library of post-processing functions
for executing a variety of wrf-python
procedures. Output is meant to be used
with the wrf-arw-tools verification
suite.

@author: aucolema
"""

import numpy as np
import xarray as xr
import wrf
from netCDF4 import Dataset
from datetime import timedelta, datetime
from wrf_ens_tools.calc import calc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from .interp_analysis import subprocess_cmd, reflectivity_to_eventgrid
import pyart
import os
from shutil import copyfile

package_dir = os.path.dirname(os.path.abspath(__file__))

dflt_var = ['td2', 'T2']
dflt_pres = [300., 500., 700., 850., 925.]
dflt_outnames = ["2m_Dewpt", "2m_Temp"]

def renameWRFOUT(enspath, ensnum, runinit,
                 subdir='mem{}/', reduced=False, ntimes=48):
    '''
    Renames wrfout files for specified ensnum. Originally
    created for use with ensemble sensitivity analysis,
    the method assumes two domains, d01 and d02, which
    become the sensitivity field domains (SENS) and
    response domains (R) respectively.

    Inputs
    ------
    enspath -- filepath string to ensemble directory
    ensnum --- int for number of ensemble members
    runinit -- datetime object for ensemble initialization time
    subdir --- string depicting member file structure in ensemble
                directory. Leave format space to format with
                member number if needed. Use empty string
                if all files are in the enspath folder, or if
                renaming for a deterministic run.
    reduced -- optional boolean specifying if wrfout files
                are reduced from their full capacity. If
                so, uses altered naming convention.
    ntimes --- optional int to specify number of forecast hours.
                Assumes 48-hr forecast unless otherwise specified.

    Outputs
    -------
    returns NULL

    Example
    -------
        enspath = '/path/to/ens/run/'
        ensnum = 42
        subdirec  = 'member{}/'
        run = datetime(2016, 5, 8, 0) ### 0 Z run on May 8, 2016
        ##### Run renameWRFOUT #####
        renameWRFOUT(enspath, ensnum, runinit=run, subdir=subdirec,
                     reduced=False, ntimes=24)

        # This file:
        # /path/to/ens/run/member1/wrfout_d01_2016-05-08_23:00:00
        # would become this file:
        # /path/to/ens/run/member1/SENS1_23.out
        # and this file:
        # /path/to/ens/run/member1/wrfout_d02_2016-05-08_23:00:00
        # would become this file:
        # /path/to/ens/run/member1/R1_23.out

        ##### If function is ran with reduced set to True #####
        renameWRFOUT(enspath, ensnum, runinit=run, subdir=subdirec,
                     reduced=True, ntimes=24)
        # Expects this file:
        # /path/to/ens/run/member1/wrfout_d02_red_2016-05-08_23:00:00
        # which becomes this file:
        # /path/to/ens/run/member1/R1_23.out
    '''
    for i in range(ensnum):
        # Format path to ensemble member of interest
        sub = enspath + subdir.format(i+1)
        for t in range(ntimes + 1):
            # Current forecast datetime object
            current = runinit + timedelta(hours = t)

            # Take care of leading zeros to build path string
            if len(str(current.month)) == 1: mo = "0" + str(current.month)
            else: mo = str(current.month)
            if len(str(current.day)) == 1: day = "0" + str(current.day)
            else: day = str(current.day)
            if len(str(current.hour)) == 1: hr = "0" + str(current.hour)
            else: hr = str(current.hour)

            # Build path string
            domain_strs = ['SENS', 'R']
            d = 1

            # If reduced, files have different naming conventions
            if reduced:
                wrfout_strs = ['{}wrfout_d0{}_red_{}-{}-{}_{}:00:00',
                               '{}wrfout_d0{}_red_{}-{}-{}_{}:00:00']
            else:
                wrfout_strs = ['{}wrfout_d0{}_{}-{}-{}_{}:00:00',
                               '{}wrfout_d0{}_{}-{}-{}_{}:00:00']

            # Assumes an outer and nested domain
            for d in range(len(domain_strs)):
                tmp = wrfout_strs[d].format(sub,
                       str(d), str(current.year), mo, day, hr)
                new = '{}{}{}_{}.out'.format(sub, domain_strs[d], str(i+1), str(t))
                tmpexists, newexists = os.path.isfile(tmp), os.path.isfile(new)
                if (tmpexists == True) and (newexists == False):
                    os.rename(tmp, new)
                d += 1
    return

# Generates dictionary of wrfout filepaths for input to process_wrf.
def gen_dict(enspath, ensnum, subdir="mem{}", ntimes=48):
    '''
    Generate dictionary of WRF outfiles sorted by
    ensemble member (or deterministic run) to feed
    into process_wrf(). Uses naming conventions
    described in renameWRFOUT().

    Inputs
    ------
    enspath -- filepath string to ensemble directory
    ensnum --- int for number of ensemble members
    subdir --- string depicting member file structure in ensemble
                directory. Leave format space to format with
                member number if needed. Use empty string
                if all files are in the enspath folder, or if
                renaming for a deterministic run.
    ntimes --- optional int to specify number of forecast hours.
                Assumes 48-hr forecast unless otherwise specified.

    Outputs
    -------
    returns dictionary of WRF outfile paths (keys are members and
    values are lists of corresponding member filepaths).
    '''
    # Initiate dictionary
    paths = {}
    for i in range(ensnum):
        # Create member index key
        key = i + 1
        # Generate list of filepaths based on assumed naming conventions
        li = ['{}{}{}'.format(enspath,
             subdir.format(i+1), '/R{}_{}.out'.format(i+1, time)) for time in range(ntimes+1)]
        # Make sure all the filepaths exist
        if np.array([os.path.exists(l) for l in li]).all():
            paths[key] = li
        else:
            raise FileNotFoundError("File(s) from member {} do not exist.".format(key))
        print(paths[key])
    return paths

# Main post-processing method for deterministic or ensemble runs.
def process_wrf(inpaths, outpath, reduced=True,
                refpath='/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREFd2',
                var=dflt_var, interp_levs=dflt_pres, outvar_names=dflt_outnames):
    '''
    A method for post-processing a list of WRF outfiles into
    one interpolated outfile for use with the verification
    suite.

    Inputs
    ------
    inpaths ------ dictionary of filepath strings for WRF outfiles
                    to process with members as keys.
    outpath ------ single filepath string for post-processed
                    WRF output.
    reduced ------ optional boolean describing whether WRF outfiles
                    are reduced from original output or complete.
                    Defaults to True, which means WRF reference
                    file will be used.
    refpath ------ optional string path to WRF reference file. Only
                    used if reduced is set to True, and defaults
                    to a reference file on Quanah that represents
                    the 4-km domain of the real-time TTU system.
    var ---------- list of strings depicting variables to process
                    and store in addition to pressure and height.
                    Defaults to dbz, slp, uh, and 10m winds variables.
                    Uses naming conventions from wrf-python module, can find at:
                    http://wrf-python.readthedocs.io/en/latest/diagnostics.html
    interp_levs -- list of floats depicting pressure levels to which heights
                    will be interpolated.

    Outputs
    -------
    returns NULL but saves post-process variables and metadata to outpath.
    '''
    # Pull original dimensions from input files
    datasets = inpaths.copy()
    for key in datasets.keys():
        datasets[key] = [Dataset(file) for file in datasets[key]]
    time = datasets[1][0].dimensions['Time']
    lat = datasets[1][0].dimensions['south_north']
    lon = datasets[1][0].dimensions['west_east']
    sigma = datasets[1][0].dimensions['bottom_top']
    # If reduced, need lats/lons and base states from reference file
    if reduced:
        wrfref = Dataset(refpath)
        wrf.disable_xarray()
        ref_cache = wrf.extract_vars(wrfref, wrf.ALL_TIMES,
                                     ("PB", "PH", "PHB", "HGT",
                                      "XLAT", "XLONG", "MAPFAC_M"))
    # Else, cache lats/lons and base state values from first file
    #  for quicker calculations.
    else:
        ref_cache = wrf.extract_vars(datasets[1][0], wrf.ALL_TIMES,
                                     ("PB", "PH", "PHB", "HGT",
                                      "XLAT", "XLONG", "MAPFAC_M"))
    # Create output file
    outfile = Dataset(outpath, 'w')
    outfile.TITLE = "OUTPUT FROM WRF-ENS-TOOLS POST PROCESS"
    outfile.START_DATE = datasets[1][0].START_DATE
    # Define dimensions
    mems = outfile.createDimension('members', len(datasets.keys()))
    outfile.createDimension(time.name, None)
    outfile.createDimension(lat.name, lat.size)
    outfile.createDimension(lon.name, lon.size)
    outfile.createDimension(sigma.name, sigma.size)

    # datasets.keys() are the member numbers, which contain
    #  all output file paths for that member in a list
    for i in range(len(datasets.keys())):
        key = list(datasets.keys())[i]
        print("Processing member/run: " + str(key))
        for t in range(len(datasets[key])):
            # Access list of filepaths for member
            files = datasets[key]
            # Process forecast hour t for member i
            dat = files[t]
            print("Processing time: " + str(t))
            # Interpolate heights, temps, and winds to pressure levels first
            wrf.disable_xarray()
            p = wrf.getvar(dat, 'p', units='hpa', cache=ref_cache)
            z = wrf.getvar(dat, 'z', units='m', cache=ref_cache)
            temp = wrf.getvar(dat, 'temp', units='degC', cache=ref_cache)
            u = wrf.getvar(dat, 'ua', units='m s-1', cache=ref_cache)
            v = wrf.getvar(dat, 'va', units='m s-1', cache=ref_cache)
            #uv = wrf.getvar(dat, 'wspd_wdir10', units='degC', cache=ref_cache)

            # Use levels from UI
            for lev in interp_levs:
                lev_ht = wrf.interplevel(z, p, lev)
                lev_temp = wrf.interplevel(temp, p, lev)
                lev_u = wrf.interplevel(u, p, lev)
                lev_v = wrf.interplevel(v, p, lev)
                ht_var = '{}_hPa_GPH'.format(int(lev))
                temp_var = '{}_hPa_T'.format(int(lev))
                uvar = '{}_hPa_U'.format(int(lev))
                vvar = '{}_hPa_V'.format(int(lev))
                if ht_var not in outfile.variables.keys():
                    outvar_hgt = outfile.createVariable(ht_var, lev_ht.dtype,
                                    (mems.name, time.name, lat.name, lon.name))
                else:
                    outvar_hgt = outfile.variables[ht_var]
                if temp_var not in outfile.variables.keys():
                    outvar_temp = outfile.createVariable(temp_var, lev_temp.dtype,
                                (mems.name, time.name, lat.name, lon.name))
                else:
                    outvar_temp = outfile.variables[temp_var]
                if uvar not in outfile.variables.keys():
                    # If uvar isn't in there, then vvar isn't either, so
                    #  create both.
                    outvar_u = outfile.createVariable(uvar, lev_u.dtype,
                                    (mems.name, time.name, lat.name, lon.name))
                    outvar_v = outfile.createVariable(vvar, lev_v.dtype,
                                    (mems.name, time.name, lat.name, lon.name))
                else:
                    outvar_u = outfile.variables[uvar]
                    outvar_v = outfile.variables[vvar]
                outvar_temp.units = 'degC'
                outvar_temp[i,t] = lev_temp[:]
                del lev_temp
                outvar_u.units = 'm s-1'
                outvar_v.units = 'm s-1'
                outvar_u[i,t] = lev_u[:]
                outvar_v[i,t] = lev_v[:]
                del lev_u, lev_v
                outvar_hgt.units = 'm'
                outvar_hgt[i,t] = lev_ht[:]
                del lev_ht
            # Process other variables requested by UI
            k = 0
            for varname in var:
                outvarnames = [outvar_names[k]]
                nvars = len(np.shape(outvarnames))
                invar = wrf.getvar(dat, varname, cache=ref_cache)
                for n in range(nvars):
                    outvarname = outvarnames[n]
                    # Create variable if not created yet
                    if outvarname not in outfile.variables.keys():
                        dimensions = [mems.name]
                        if reduced:
                            # Enable xarray to get units of the variable
                            wrf.enable_xarray()
                            dim_var = wrf.getvar(wrfref, varname)
                            units = dim_var.units
                            # Create list of dimensions
                            for d in dim_var.dims:
                                dimensions.append(d)
                                if d not in outfile.dimensions:
                                    outfile.createDimension(d, dim_var.sizes[d])
                            # If Time not the second dimension, insert it
                            if dimensions[1] != t:
                                dimensions.insert(1, time.name)
                            wrf.disable_xarray()
                        else:
                            # Xarray already enabled, so pull units directly
                            units = invar.units
                            # Create list of dimensions
                            for d in invar.dims:
                                dimensions.append(d)
                                if d not in outfile.dimensions:
                                    outfile.createDimension(d, invar.sizes[d])
                            # If Time not the second dimension, insert it
                            if dimensions[1] != t:
                                dimensions.insert(1, time.name)
                        # Create variable and assign units
                        tuple(dimensions)
                        outvar = outfile.createVariable(outvarname, invar.dtype, dimensions)
                        outvar.units = units
                    # Otherwise, pull variable from outfile
                    else:
                        outvar = outfile.variables[outvarname]
                    if len(np.shape(invar)) > 2:
                        outvar[i, t] = invar[n]
                    else:
                        outvar[i, t] = invar[:]
                    del invar
                k += 1
            dat.close()
        outfile.close()
        return

def postTTUWRFanalysis(inpath, outpath,
                       refpath='/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREF'):
    """
    Post-process Brian's TTU WRF analyses into wrf-ens-tools
    naming conventions and store to new netCDF to use in
    subsetting.

    Here's how Brian stratifies the variables:
      t(:,:,1)=gph300next(:,:)
      t(:,:,2)=gph500next(:,:)
      t(:,:,3)=gph700next(:,:)
      t(:,:,4)=gph850next(:,:)
      t(:,:,5)=gph925next(:,:)
      t(:,:,6)=t300next(:,:)
      t(:,:,7)=t500next(:,:)
      t(:,:,8)=t700next(:,:)
      t(:,:,9)=t850next(:,:)
      t(:,:,10)=t925next(:,:)
      t(:,:,11)=u300next(:,:)
      t(:,:,12)=u500next(:,:)
      t(:,:,13)=u700next(:,:)
      t(:,:,14)=u850next(:,:)
      t(:,:,15)=u925next(:,:)
      t(:,:,16)=v300next(:,:)
      t(:,:,17)=v500next(:,:)
      t(:,:,18)=v700next(:,:)
      t(:,:,19)=v850next(:,:)
      t(:,:,20)=v925next(:,:)
      t(:,:,21)=td300next(:,:)
      t(:,:,22)=td500next(:,:)
      t(:,:,23)=td700next(:,:)
      t(:,:,24)=td850next(:,:)
      t(:,:,25)=td925next(:,:)
      t(:,:,26)=q300next(:,:)
      t(:,:,27)=q500next(:,:)
      t(:,:,28)=q700next(:,:)
      t(:,:,29)=q850next(:,:)
      t(:,:,30)=q925next(:,:)
      t(:,:,31)=slpnext(:,:)
      t(:,:,32)=t2next(:,:)
      t(:,:,33)=td2next(:,:)
      t(:,:,34)=u10next(:,:)
      t(:,:,35)=v10next(:,:)
      t(:,:,36)=0.0
      t(:,:,37)=0.0

    Convert to following naming conventions in a new netCDF:
    sensstringslist = ["300 hPa GPH","500 hPa GPH","700 hPa GPH",
                       "850 hPa GPH","300 hPa T","500 hPa T",
                       "700 hPa T","850 hPa T","925 hPa T",
                       "300 hPa U-Wind","500 hPa U-Wind",
                       "700 hPa U-Wind","850 hPa U-Wind",
                       "925 hPa U-Wind","300 hPa V-Wind","500 hPa V-Wind",
                       "700 hPa V-Wind","850 hPa V-Wind","925 hPa V-Wind",
                       "SLP","2m Temp","2m Q",
                       "2m Dewpt","10m U-Wind","10m V-Wind"]
    """
    # Pull analysis variables
    og_analysis = Dataset(inpath)
    anlvars = og_analysis.variables['T'][0]
    gph300 = anlvars[0,:,:]
    gph500 = anlvars[1,:,:]
    gph700 = anlvars[2,:,:]
    gph850 = anlvars[3,:,:]
    gph925 = anlvars[4,:,:]
    temp300 = anlvars[5,:,:]
    temp500 = anlvars[6,:,:]
    temp700 = anlvars[7,:,:]
    temp850 = anlvars[8,:,:]
    temp925 = anlvars[9,:,:]
    u300 = anlvars[10,:,:]
    u500 = anlvars[11,:,:]
    u700 = anlvars[12,:,:]
    u850 = anlvars[13,:,:]
    u925 = anlvars[14,:,:]
    v300 = anlvars[15,:,:]
    v500 = anlvars[16,:,:]
    v700 = anlvars[17,:,:]
    v850 = anlvars[18,:,:]
    v925 = anlvars[19,:,:]
    td300 = anlvars[20,:,:]
    td500 = anlvars[21,:,:]
    td700 = anlvars[22,:,:]
    td850 = anlvars[23,:,:]
    td925 = anlvars[24,:,:]
    q300 = anlvars[25,:,:]
    q500 = anlvars[26,:,:]
    q700 = anlvars[27,:,:]
    q850 = anlvars[28,:,:]
    q925 = anlvars[29,:,:]
    slp = anlvars[30,:,:]
    t2 = anlvars[31,:,:]
    td2 = anlvars[32,:,:]
    u10 = anlvars[33,:,:]
    v10 = anlvars[34,:,:]

    wrf_d1 = Dataset(refpath)
    lons, lats = wrf_d1.variables['XLONG'][0], wrf_d1.variables['XLAT'][0]
    wrf_idim = len(lons[0,:])
    wrf_jdim = len(lats[:,0])

    sensvarlist = [gph300,gph500,gph700,gph850,gph925,temp300,temp500,temp700,
                   temp850,temp925,u300,u500,u700,u850,u925,v300,
                   v500,v700,v850,v925,td300,td500,td700,td850,
                   td925,q300,q500,q700,q850,q925,slp,t2,td2,u10,v10]
    sensstringslist = ["300 hPa GPH","500 hPa GPH","700 hPa GPH",
                       "850 hPa GPH","925 hPa GPH","300 hPa T","500 hPa T",
                       "700 hPa T","850 hPa T","925 hPa T","300 hPa U-Wind",
                       "500 hPa U-Wind","700 hPa U-Wind","850 hPa U-Wind",
                       "925 hPa U-Wind","300 hPa V-Wind","500 hPa V-Wind",
                       "700 hPa V-Wind","850 hPa V-Wind","925 hPa V-Wind",
                       "300 hPa Dewpt", "500 hPa Dewpt", "700 hPa Dewpt",
                       "850 hPa Dewpt", "925 hPa Dewpt", "300 hPa Q",
                       "500 hPa Q", "700 hPa Q", "850 hPa Q", "925 hPa Q",
                       "SLP","2m Temp","2m Dewpt",
                       "10m U-Wind","10m V-Wind"]


    # Write interpolated variables to netCDF
    new_analysis = Dataset(outpath, "w", format="NETCDF4")
    new_analysis.createDimension('lat', wrf_jdim)
    new_analysis.createDimension('lon', wrf_idim)
    new_analysis.createDimension('time', None)
    xlat = new_analysis.createVariable("XLAT", float, dimensions=('lat','lon'))
    xlat[:,:] = lats
    xlon = new_analysis.createVariable("XLONG", float, dimensions=('lat','lon'))
    xlon[:,:] = lons

    # Interpolate and save!!
    for i in range(len(sensvarlist)):
        var = new_analysis.createVariable(sensstringslist[i].replace(" ","_"),
                                          sensvarlist[i].dtype,
                                          dimensions=('lat','lon'))
        var[:,:] = sensvarlist[i]
    new_analysis.close()

# Post-process practically perfect
def storePracPerfNativeGrid(modelinit, fcsthrs, outpath, nbrhd, dx, sixhour=False):
    '''
    Calculate hourly or six-hourly practically perfect on WRF grid.
    Save to netCDF file for verification.

    Inputs
    ------
    modelinit - WRF run initialization datetime obj
                for formatting times correctly.
                (WRF run doesn't need to exist,
                just need baseline datetime to add
                forecast hours to.)
    fsthrs ---- list of forecast hour integers, or
                hours since modelinit time.
    sixhour --- boolean describing whether to store
                practically perfect probs in six-hour or
                one-hour increments. Defaults to one hour.
    outpath --- string specifying absolute path of
                netCDF output.
    nbrhd ----- neighborhoood distance in km to use for
                the distance-based sigma of the Gaussian
                kernel.
    dx -------- horizontal grid-spacing on model domain to
                serve as the denominator of the Gaussian
                kernel.

    Outputs
    -------
    returns NULL, but saves to netCDF outpath.
    '''
    # Create outfile
    netcdf_out = Dataset(outpath, "w", format="NETCDF4")
    # Calculate pperf to pull lat/lon data
    sigma = nbrhd / dx
    pperf, lon, lat = calc_prac_perf_native_grid(modelinit, sixhour,
                                     fcsthrs[0], sigma=sigma)
    # Set up netCDF
    netcdf_out.START_DATE = modelinit.strftime('%Y-%m-%d_%H:%M:%S')
    netcdf_out.createDimension('Time', len(fcsthrs))
    netcdf_out.createDimension('south_north', len(lat[:,0]))
    netcdf_out.createDimension('west_east', len(lon[0,:]))
    netcdf_out.createDimension('sigma', 1)
    times = netcdf_out.createVariable('fhr', int, ('Time'))
    sig = netcdf_out.createVariable('sigma', float, ('sigma'))
    pperfout = netcdf_out.createVariable('practically_perfect', float, ('Time', 'south_north', 'west_east'))
    sig[:] = sigma
    # Populate outfile with pperf
    for t in range(len(fcsthrs)):
        pperf, lon, lat = calc_prac_perf_native_grid(modelinit,
                                sixhour, fcsthrs[t], sigma=sigma)
        pperfout[t] = pperf[:]*100.
        times[t] = fcsthrs[t]
    netcdf_out.close()
    return

def storePracPerfSPCGrid(modelinit, fcsthrs, outpath,
                            nbrhd, dx, sixhour=False,
                            wrfrefpath='/lustre/scratch/aucolema/2016052600/wrfoutREFd2'):
    '''
    Calculate hourly or six-hourly practically perfect on 80-km grid-spacing
    SPC grid, which is then interpolated to WRF grid and stored to netCDF
    for verification.

    Inputs
    ------
    modelinit - WRF run initialization datetime obj
                for formatting times correctly.
                (WRF run doesn't need to exist,
                just need baseline datetime to add
                forecast hours to.)
    fsthrs ---- list of forecast hour integers, or
                hours since modelinit time.
    sixhour --- boolean describing whether to store
                practically perfect probs in six-hour or
                one-hour increments. Defaults to one hour.
    outpath --- string specifying absolute path of
                netCDF output.
    nbrhd ----- neighborhoood distance in km to use for
                the distance-based sigma of the Gaussian
                kernel.
    dx -------- horizontal grid-spacing on model domain to
                serve as the denominator of the Gaussian
                kernel.

    Outputs
    -------
    returns NULL, but saves to netCDF outpath.
    '''
    # Create outfile
    netcdf_out = Dataset(outpath, "w", format="NETCDF4")
    # Calculate pperf to pull lat/lon data
    sigma = nbrhd / dx
    pperf, lon, lat = calc.calc_prac_perf_spc_grid(modelinit, sixhour,
                                     fcsthrs[0], sigma=sigma,
                                     wrfrefpath=wrfrefpath)
    # Set up netCDF
    netcdf_out.START_DATE = modelinit.strftime('%Y-%m-%d_%H:%M:%S')
    netcdf_out.createDimension('Time', len(fcsthrs))
    netcdf_out.createDimension('south_north', len(lat[:,0]))
    netcdf_out.createDimension('west_east', len(lon[0,:]))
    netcdf_out.createDimension('sigma', 1)
    times = netcdf_out.createVariable('fhr', int, ('Time'))
    sig = netcdf_out.createVariable('sigma', float, ('sigma'))
    pperfout = netcdf_out.createVariable('practically_perfect', float, ('Time', 'south_north', 'west_east'))
    sig[:] = sigma
    # Populate outfile with pperf
    for t in range(len(fcsthrs)):
        pperf, lon, lat = calc.calc_prac_perf_spc_grid(modelinit, sixhour, fcsthrs[t], sigma=sigma)
        pperfout[t] = pperf[:]*100.
        times[t] = fcsthrs[t]
    netcdf_out.close()
    return

# Interpolate storm reports to nearest grid point - mainly for reliability calc
def storeNearestNeighbor(modelinit, fcsthrs, outpath, sixhour=True,
                            wrfrefpath='/lustre/research/bancell/aucolema/HWT2016runs'
                                       '/2016050800/wrfoutREFd2'):
    '''
    Calculate hourly nearest neighbor on WRF grid.
    Save to netCDF file for verification.

    Inputs
    ------
    modelinit - WRF run initialization datetime obj
                for formatting times correctly.
                (WRF run doesn't need to exist,
                just need baseline datetime to add
                forecast hours to.)
    fsthrs ---- list of forecast hour integers, or
                hours since modelinit time.
    outpath --- string specifying absolute path of
                netCDF output.
    sixhour --- boolean specifying whether to store
                six-hour or one-hour period of
                nearest neighbor storm reports

    Outputs
    -------
    returns NULL, but saves to netCDF outpath.
    '''
    # Create outfile
    netcdf_out = Dataset(outpath, "w", format="NETCDF4")
    # Calculate pperf to pull lat/lon data
    grid = calc.nearest_neighbor_spc(modelinit, sixhour, fcsthrs[0], nbrhd=0.,
                                wrfrefpath=wrfrefpath)
    # Set up netCDF
    netcdf_out.START_DATE = modelinit.strftime('%Y-%m-%d_%H:%M:%S')
    netcdf_out.SIXHOUR = str(sixhour)
    netcdf_out.createDimension('Time', len(fcsthrs))
    netcdf_out.createDimension('south_north', len(grid[:,0]))
    netcdf_out.createDimension('west_east', len(grid[0,:]))
    times = netcdf_out.createVariable('fhr', int, ('Time'))
    nearest_out = netcdf_out.createVariable('nearest_neighbor', float, ('Time',
                                            'south_north', 'west_east'))
    # Populate outfile with pperf
    for t in range(len(fcsthrs)):
        grid = calc.nearest_neighbor_spc(modelinit, sixhour, fcsthrs[t], nbrhd=0.,
                                    wrfrefpath=wrfrefpath)
        nearest_out[t] = grid[:]
        times[t] = fcsthrs[t]
    netcdf_out.close()
    return

def storeNearestNeighborFortran(modelinit, fcsthr, outpath, sixhour=True,
                            variable="updraft_helicity",
                            wrfrefpath='/lustre/research/bancell/aucolema/HWT2016runs'
                                       '/2016050800/wrfoutREFd2',
                            interpgridradfiles=None,
                            reflthreshold=None):
    '''
    Calculate hourly nearest neighbor on WRF grid.
    Overwrite to WRF outfile for fast verification
    with Fortran 77.

    Inputs
    ------
    modelinit ------------- WRF run initialization datetime obj
                            for formatting times correctly.
                            (WRF run doesn't need to exist,
                            just need baseline datetime to add
                            forecast hours to.)
    fcsthr ---------------- forecast hour integer in
                            hours since modelinit time.
    outpath --------------- string specifying absolute path of
                            netCDF output.
    sixhour --------------- boolean specifying whether to store
                            six-hour or one-hour period of
                            nearest neighbor storm reports.
    variable -------------- defaults to 'updraft_helicity' but
                            also supports reliability processing
                            for 'reflectivity'.
    interpgridradfiles ---- if verifying reflectivity, will need
                            to provide list of interpolated
                            GridRad radar filepaths.
    reflthreshold --------- if verifying reflectivity, will need
                            to provide an exceedance threshold as
                            a float.

    Outputs
    -------
    returns NULL, but saves to netCDF outpath.
    '''
    # Create outfile
    #os.popen("cp {} {}".format(wrfrefpath, outpath))
    # os.popen("cp {} {}".format(wrfrefpath, outpath))
    copyfile(wrfrefpath, outpath)
    netcdf_out = Dataset(outpath, "a")
    obvar = "P_HYD"
    netcdf_out.START_DATE = modelinit.strftime('%Y-%m-%d_%H:%M:%S')
    netcdf_out.SIXHOUR = str(sixhour)
    netcdf_out.FHR = fcsthr

    # Populate outfile with gridded observations
    if variable == "updraft_helicity":
        grid = calc.nearest_neighbor_spc(modelinit, sixhour, fcsthr, nbrhd=0.,
                                    wrfrefpath=wrfrefpath)
    elif variable == "reflectivity":
        grid = reflectivity_to_eventgrid(interp_gridradfiles=interpgridradfiles,
                runinitdate=modelinit, sixhr=sixhour, rtime=fcsthr,
                threshold=reflthreshold)
    netcdf_out.variables[obvar][0,0,:,:] = grid[:,:]

    netcdf_out.close()
    return

def storeReliabilityRboxFortran(basedir, fcsthr, probpath, obpath, outfile,
                rboxpath, sixhour=True, variable="updraft_helicity",
                rthresh=25.0, nbrhd=30.0,
                wrfrefpath='/lustre/research/bancell/aucolema/HWT2016runs'
                '/2016050800/wrfoutREFd2'):
    """
    Calculates reliability using fortran 77 code and stores
    it to 'outfile' in the same directory in which the function was called.

    Inputs
    ------
    basedir ----------- absolute path to the base directory
                        in which the reliability output file
                        will be stored.
    fcsthr ------------ integer in number of forecast hours
                        since run initialization.
    probpath ---------- absolute filepath to probability
                        netCDF file.
    obpath ------------ absolute path to output file from
                        storeNearestNeighborFortran()
                        call.
    outfile ----------- absolute path to file where
                        reliability will be stored.
    rboxpath ---------- absolute path to "esens.in" file
                        to pull response box bounds from.
    nbrhd ------------- searching radius in km to calculate
                        reliability with (should be the same
                        neighborhood used in probability
                        calculations).

    Outputs
    -------
    returns NULL but stores all response box reliability
    statistics to 'reliability_out.nc' in the given
    base directory.
    """
    esensin = np.genfromtxt(rboxpath)
    rbox_bounds = esensin[4:8]
    if os.path.exists(basedir+'reliability.in'):
        subprocess_cmd("rm {}/reliability.in".format(basedir))
    args = [str("\\'"+probpath+"\\'"), str("\\'"+obpath+"\\'"),
            str("\\'"+outfile+"\\'"),
            fcsthr, str("\\'"+variable+"\\'"),
            rthresh, sixhour, nbrhd, rbox_bounds[0],
            rbox_bounds[1], rbox_bounds[2], rbox_bounds[3]]
    print("Throwing these into reliability.in")
    print(args)
    # np.savetxt(basedir+"reliability.in", np.asarray(args))
    os.popen("echo {} > {}reliability.in".format(args[0], basedir))
    for arg in args[1::]:
        os.popen("echo {} >> {}reliability.in".format(arg, basedir))
    if os.path.exists(outfile):
        subprocess_cmd("rm {}".format(outfile))
        print("Removed old reliability output file...")
    print("Running reliability arguments...")
    print("Directory:", package_dir)
    subprocess_cmd("{}/reliabilitycalc <{}/reliability.in \
                    >{}reliability.out".format(package_dir, basedir, basedir))
    return

# TO-DO: to make real-time useable, add prob calculation
# with procalcSUBSET.f
def calcEnsProbs(enspath, members, wrfref):
    pass

################### GridRad Reflectivity post-processing #######################
def merge_refl_data(gridradfiles):
    """
    Merges given list of GridRad reflectivity data over the time dimension
    and returns a concatenated xarray DataArray object.
    """
    # Read data
    dats = [xr.open_dataset(file) for file in gridradfiles]
    inds = dats[0]["index"]
    time = dats[0].Analysis_time
    lats = dats[0]["Latitude"]
    lons = dats[0]["Longitude"]-360.
    z = dats[0]["Altitude"]
    refl = dats[0]["Reflectivity"].values
    # Reshape reflectivity data
    refl_vals = np.zeros(len(z.values)*len(lats.values)*len(lons.values))*np.NaN
    refl_vals[inds.values] = refl
    refl_reshape = refl_vals.reshape((len(z.values), len(lats.values),
                                        len(lons.values)))
    # Create reflectivity xarray DataArray obj
    dat = xr.DataArray(refl_reshape, coords=[z, lats, lons])
    dat.name = "Reflectivity"
    # For each time stamp in the dataset, pull reflectivity and concatenate
    #  onto DataArray object
    for i in range(1,len(dats)):
        # Pull reflectivity information and post-process
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
        # Initialize DataArray to concatenate
        da = xr.DataArray(refl_reshape, coords=[z, lats, lons])
        da.attrs['Analysis_Endtime'] = time
        da.name = "Reflectivity"
        # Concatenate new DataArray onto original
        dat = xr.concat([dat, da], dim='Hours')
    # Add the end datetime object as an attribute
    dat.attrs["Analysis_Endtime"] = time
    return dat

# Gridrad subset data to a response box
def plot_refl_rbox(gridradfiles, rboxpath, zlev):
    """
    Read a GridRad file over a specified portion of the
    country, as specified by an esens.in file.

    Inputs
    ------
    infile -------- absolute path to GridRad nc file to
                    process
    rboxpath ------ absolute path to esens.in file provided
                    to sensitivity code that contains
                    response box bounds
    zlev ---------- zero-based vertical level to plot

    Outputs
    -------
    returns data dictionary similar to read_file() output
    except only over a response box
    """
    # Read response box bounds
    esensin = np.genfromtxt(rboxpath)
    rbox_bounds = esensin[4:8]

    # Merge data
    dat = merge_refl_data(gridradfiles)
    lons = dat["Longitude"]
    lats = dat["Latitude"]

    # Build meshgrid
    xmesh, ymesh = np.meshgrid(lons, lats)

    ########### Plot original radar data over response box ##################
    # Build response box
    llon, ulon, llat, ulat = rbox_bounds
    width = ulon - llon
    height = ulat - llat

    # Initialize plot
    fig = plt.figure(figsize=(10, 10))

    # Build projection/map
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
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

    for i in range(len(dat.values)):
        endtime_minus_nhrs = len(dat.values) - i
        time = str(datetime.strptime(dat.Analysis_Endtime, "%Y-%m-%d %H:%M:%SZ") + \
                timedelta(hours=-1*endtime_minus_nhrs)).replace(" ", "_")
        print("Plotting GridRad data for", time)
        # Plot radar data
        refl = ax.contourf(xmesh, ymesh, dat.values[i,zlev],
                    transform=ccrs.PlateCarree(), cmap="pyart_HomeyerRainbow")
        plt.colorbar(refl, ax=ax, label="Reflectivity",
                    fraction=0.0289, pad=0.0)
        plt.title("GridRad Vertical-Level-{} Reflectivity Data valid {}".format(zlev+1,
                    time))
        figname = 'gridrad_zlev{}_valid{}.png'.format(zlev, time)
        plt.savefig(figname)
    return
