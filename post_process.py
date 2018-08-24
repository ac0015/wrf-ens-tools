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
import wrf
from netCDF4 import Dataset
from datetime import timedelta
from calc import calc_prac_perf
import os

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
                #print(outvarnames)
                nvars = len(np.shape(outvarnames))
                #print(varname, dat)
                invar = wrf.getvar(dat, varname, cache=ref_cache)
                for n in range(nvars): 
                    outvarname = outvarnames[n]
                    #print(outvarname)
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

# Post-process practically perfect
def storePracPerf(modelinit, fcsthrs, outpath, sigma=2):
    '''
    Calculate hourly practically perfect on WRF grid.
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

    Outputs
    -------
    returns NULL, but saves to netCDF outpath. 
    '''
    # Create outfile
    netcdf_out = Dataset(outpath, 'w')
    # Calculate pperf to pull lat/lon data
    pperf, lon, lat = calc_prac_perf(modelinit, False, fcsthrs[0], sigma=sigma)
    # Set up netCDF
    netcdf_out.START_DATE = modelinit.strftime('%Y-%m-%d_%H:%M:%S')
    netcdf_out.createDimension('Time', len(fcsthrs))
    netcdf_out.createDimension('south_north', len(lat[:,0]))
    netcdf_out.createDimension('west_east', len(lon[0,:]))
    netcdf_out.createDimension('sigma', 1)
    times = netcdf_out.createVariable('fhr', int, ('Time'))
    sig = netcdf_out.createVariable('sigma', int, ('sigma'))
    pperfout = netcdf_out.createVariable('practically_perfect', float, ('Time', 'south_north', 'west_east'))
    sig[:] = sigma
    # Populate outfile with pperf
    for t in range(len(fcsthrs)):
        pperf, lon, lat = calc_prac_perf(modelinit, False, fcsthrs[t], sigma=sigma)
        pperfout[t] = pperf[:]*100.
        times[t] = fcsthrs[t]
    netcdf_out.close()
    return
     
# TO-DO: to make real-time useable, add prob calculation
# with procalcSUBSET.f
def calcEnsProbs(enspath, members, wrfref):
    pass
                