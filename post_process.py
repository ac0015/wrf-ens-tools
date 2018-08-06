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
from datetime import datetime, timedelta
from shutil import copyfile
from calc import calc_prac_perf
import os

dflt_var = ['dbz', 'slp', 'updraft_helicity', 'wspd_wdir10']
dflt_pres = [300., 500., 700., 850., 925.]

def gen_dict(enspath, nmems, ntimes, subdir="mem{}"):
    '''
    Generate dictionary of WRF outfiles sorted by
    ensemble member (or deterministic run) to feed 
    into process_wrf(). Uses Austin's naming conventions.
    '''
    paths = {}
    for i in range(nmems):
        key = i + 1
        li = ['{}{}{}'.format(enspath, 
             subdir.format(i+1), '/R{}_{}.out'.format(i+1, time)) for time in range(ntimes+1)]
        if np.array([os.path.exists(l) for l in li]).all():
            paths[key] = li
        else:
            raise FileNotFoundError("File(s) from member {} do not exist.".format(key))    
        print(paths[key])
    return paths

# Main post-processing method for deterministic or ensemble runs.
def process_wrf(inpaths, outpath, reduced=True, 
                refpath='/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREFd2',
                var=dflt_var, interp_levs=dflt_pres):
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
    
    '''
    # Pull original dimensions from input files
    datasets = inpaths.copy()
    for key in datasets.keys():
        datasets[key] = [Dataset(file) for file in datasets[key]]
    time = datasets[1][0].dimensions['Time']
    lat = datasets[1][0].dimensions['south_north']
    lon = datasets[1][0].dimensions['west_east']
    sigma = datasets[1][0].dimensions['bottom_top']
    # If reduced, need lats and lons from reference file
    if reduced:
        wrfref = Dataset(refpath)
        wrf.disable_xarray()
        ref_cache = wrf.extract_vars(wrfref, wrf.ALL_TIMES,  
                                     ("PB", "PH", "PHB", "HGT", 
                                      "XLAT", "XLONG", "MAPFAC_M"))
    else:
        ref_cache = wrf.extract_vars(datasets[1][0], wrf.ALL_TIMES,  
                                     ("PB", "PH", "PHB", "HGT", 
                                      "XLAT", "XLONG", "MAPFAC_M"))
    # Create output file
    outfile = Dataset(outpath, 'w')
    outfile.TITLE = "OUTPUT FROM WRF-ARW-TOOLS POST PROCESS"
    outfile.START_DATE = datasets[1][0].START_DATE
    # Define dimensions
    mems = outfile.createDimension('members', len(datasets.keys()))
    outfile.createDimension(time.name, None)
    outfile.createDimension(lat.name, lat.size)
    outfile.createDimension(lon.name, lon.size)
    outfile.createDimension(sigma.name, sigma.size)
    
    for i in range(len(datasets.keys())):
        key = list(datasets.keys())[i]
        print("Processing member/run: " + str(key))
        for t in range(len(datasets[key])):
            files = datasets[key]
            dat = files[t]
            print("Processing time: " + str(t))
            # Interpolate heights to pressure levels first
            p = wrf.getvar(dat, 'p', units='hpa', cache=ref_cache)
            z = wrf.getvar(dat, 'z', units='m', cache=ref_cache)
            # Use levels from UI
            for lev in interp_levs:
                lev_ht = wrf.interplevel(z, p, lev)
                ht_var = 'hgt_{}'.format(int(lev))
                if ht_var not in outfile.variables.keys():
                    outvar = outfile.createVariable(ht_var, lev_ht.dtype, 
                                    (mems.name, time.name, lat.name, lon.name))
                else:
                    outvar = outfile.variables[ht_var]
                outvar.units = 'm'
                outvar[i,t] = lev_ht[:]
                del lev_ht
            for varname in var:
                invar = wrf.getvar(dat, varname, cache=ref_cache)
                if varname not in outfile.variables.keys():
                    dimensions = [mems.name]
                    if reduced:
                        wrf.enable_xarray()
                        dim_var = wrf.getvar(wrfref, varname)
                        units = dim_var.units
                        for d in dim_var.dims:
                            dimensions.append(d)
                            if d not in outfile.dimensions:
                                outfile.createDimension(d, dim_var.sizes[d])
                        if dimensions[1] != t:
                            dimensions.insert(1, time.name)
                        wrf.disable_xarray()
                    else:
                        units = invar.units
                        for d in invar.dims:
                            dimensions.append(d)
                            if d not in outfile.dimensions:
                                outfile.createDimension(d, invar.sizes[d])
                        if dimensions[1] != t:
                            dimensions.insert(1, time.name)
                    tuple(dimensions)
                    outvar = outfile.createVariable(varname, invar.dtype, dimensions)
                    outvar.units = units
                else:
                    outvar = outfile.variables[varname]
                outvar[i, t] = invar[:]
                del invar
            dat.close()
    outfile.close()
    return

# Post-process practically perfect
def storePracPerf(modelinit, fcsthrs, outpath):
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
    Returns nothing, but saves to netCDF outpath. 
    '''
    # Create outfile
    netcdf_out = Dataset(outpath, 'w')
    # Calculate pperf to pull lat/lon data
    pperf, lon, lat = calc_prac_perf(modelinit, False, fcsthrs[0])
    # Set up netCDF
    netcdf_out.START_DATE = modelinit.strftime('%Y-%m-%d_%H:%M:%S')
    netcdf_out.createDimension('Time', len(fcsthrs))
    netcdf_out.createDimension('south_north', len(lat[:,0]))
    netcdf_out.createDimension('west_east', len(lon[0,:]))
    times = netcdf_out.createVariable('fhr', int, ('Time'))
    pperfout = netcdf_out.createVariable('practically_perfect', float, ('Time', 'south_north', 'west_east'))
    # Populate outfile with pperf
    for t in range(len(fcsthrs)):
        pperf, lon, lat = calc_prac_perf(modelinit, False, fcsthrs[t])
        pperfout[t] = pperf[:]*100.
        times[t] = fcsthrs[t]
    netcdf_out.close()
    return
     
# TO-DO: to make real-time useable, add prob calculation
# with procalcSUBSET.f
def calcEnsProbs(enspath, members, wrfref):
    pass
                