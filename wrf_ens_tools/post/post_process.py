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
from sklearn.neighbors import BallTree
import dask.array as da
import xarray as xr
import metpy.constants as constants
import wrf
from netCDF4 import Dataset
from datetime import timedelta, datetime
from wrf_ens_tools.calc import calc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pyproj
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from .interp_analysis import subprocess_cmd, reflectivity_to_eventgrid, bilinear_interp
# import pyart
import os
from shutil import copyfile
import warnings

package_dir = os.path.dirname(os.path.abspath(__file__))

P0 = constants.P0.to('Pa').m
kappa = constants.kappa.m

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
    return

def postIdealizedAnalysis(inpath, outpath, member,
                        refpath='/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREF'):
    """
    Siphons post-processed sensitivity variable values from a "SENSvals.nc" or
    similar file produced by the sensvector executable found in the sensitivity
    module. Stores the specified member's post-processed fields in a specified
    outpath. Intended to be used in idealized testing experiments where one
    ensemble member is randomly-selected to serve as truth.

    Inputs
    ------
    inpath ------------ input path to netCDF file where post-processed member
                        fields are stored (valid at the sensitivity time).
                        Should be formatted as output from sensvector exec.
    outpath ----------- output path where member sensitivity variable fields
                        will be stored like an analysis.
    member ------------ one-based integer specifying the ensemble member to
                        pull from the input file.
    refpath ----------- absolute path to reference file containing geographical
                        data and other relevant metadata.

    Outputs
    -------
    Returns NULL, but stores member sensitivity variable data valid at the
    sensitivity time in the specified outpath.
    """
    # SENSvals file naming conventions
    sensval_varstrings = ["GPH_300", "GPH_500", "GPH_700", "GPH_850", "SKIP",
                            "T_300", "T_500", "T_700", "T_850", "T_925",
                            "U_300", "U_500", "U_700", "U_850", "U_925",
                            "V_300", "V_500", "V_700", "V_850", "V_925",
                            "SKIP", "SKIP", "SKIP", "SKIP", "SKIP", "SKIP",
                            "SKIP", "SKIP", "Q_850", "SKIP", "SLP", "T2",
                            "TD2", "U10", "V10"]
    # Post-processed new file naming conventions
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

    # Get more dimensions/geographical info
    wrf_d1 = Dataset(refpath)
    lons, lats = wrf_d1.variables['XLONG'][0], wrf_d1.variables['XLAT'][0]
    wrf_idim = len(lons[0,:])
    wrf_jdim = len(lats[:,0])

    # Write interpolated variables to netCDF
    new_analysis = Dataset(outpath, "w", format="NETCDF4")
    new_analysis.createDimension('lat', wrf_jdim)
    new_analysis.createDimension('lon', wrf_idim)
    new_analysis.createDimension('time', None)
    xlat = new_analysis.createVariable("XLAT", float, dimensions=('lat','lon'))
    xlat[:,:] = lats
    xlon = new_analysis.createVariable("XLONG", float, dimensions=('lat','lon'))
    xlon[:,:] = lons

    # Open dataset and start pulling member fields
    member_fields = np.zeros((len(sensval_varstrings), wrf_jdim, wrf_idim))
    sensvar_dat = Dataset(inpath)
    for ind, var in enumerate(sensval_varstrings):
        # print("SENSvals variable:", var, "New variable string", sensstringslist[ind])
        if var != "SKIP":
            member_fields[ind] = sensvar_dat[var][member-1][:]
            newvar = new_analysis.createVariable(
                                        sensstringslist[ind].replace(" ","_"),
                                        member_fields[ind].dtype,
                                        dimensions=('lat','lon'))
            newvar[:,:] = member_fields[ind]
    new_analysis.close()
    return

# Post-process practically perfect
def storePracPerfNativeGrid(modelinit, fcsthrs, outpath, nbrhd, dx,
                            sixhour=False):
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
    Calculate hourly or six-hourly practically perfect on 81-km grid-spacing
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
    idealized - optional boolean describing whether to
                accept an SSPF array-like to store to
                practically-perfect file or to calculate
                practically-perfect probs with actual
                storm reports from the specified date.
    SSPF ------ if idealized is set to True, will need to
                provide array-like with surrogate severe
                probabilities to used as idealized PP probs.


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
    # Populate outfile with pperf (if six-hour, calculating six-hr pperf for each fcst hr)
    for t in range(len(fcsthrs)):
        pperf, lon, lat = calc.calc_prac_perf_spc_grid(modelinit, sixhour, fcsthrs[t], sigma=sigma)
        pperfout[t] = pperf[:]*100.
        times[t] = fcsthrs[t]
    netcdf_out.close()
    return

def storeIdealizedPracPef(sspf_arr, outlats, outlons, outpath,
                            sigma, modelinit, fhrs, spc_grid=True):
    """
    Takes array of pre-calculated surrogate severe probability forecast
    and meta data pertaining to that SSPF and stores it to a specified
    netCDF outpath. If spc_grid is set to True, will interpolate from
    lat/lons specified in pperf_grid_template.npz to lat/lons provided
    as input before storing to output file.

    Inputs
    ------
    sspf_arr ---------- array-like containing surrogate severe probability
                        forecast with values from 0-1
    outlats ----------- 2D array-like containing desired output latitudes
                        for which sspf_arr will be valid. If spc_grid=False,
                        this lat array should correspond to the sspf_arr. If
                        spc_grid=True, this lat array should correspond with
                        the latitudes that the sspf_arr will be interpolated
                        to
    outlons ----------- 2D array-like containing desired output longitudes
                        for which sspf_arr will be valid. If spc_grid=False,
                        this longitude array should correspond with the
                        longitudes that the sspf_arr will be interpolated to
    outpath ----------- absolute filepath to describing location to store
                        practically-perfect output
    sigma ------------- smoothing parameter (nbrhd/dx) with which
                        Guassian kernel was applied (for meta data)
    modelinit --------- datetime object dictating the run initialization
                        used
    fhrs -------------- array-like of forecast hours for which pperf is
                        valid
    spc_grid ---------- optional boolean describing whether original data
                        is stored on the SPC 211 grid, in which case data
                        will be interpolated to outlats/outlons

    Output
    ------
    returns NULL but stores PP probs to specified output filepath
    """
    # Create outfile
    netcdf_out = Dataset(outpath, "w", format="NETCDF4")

    # If spc_grid=True, pull lat/lons from pperf grid and interpolate
    if spc_grid:
        print("Interpolating practically-perfect to SPC grid")
        f = np.load(package_dir + "/pperf_grid_template.npz")
        lats = f['lat']
        lons = f['lon']
        sspf_arr = bilinear_interp(grid1x=lons, grid1y=lats,
                        grid2x=outlons, grid2y=outlats, z=sspf_arr)

    # Set up netCDF
    netcdf_out.START_DATE = modelinit.strftime('%Y-%m-%d_%H:%M:%S')
    netcdf_out.createDimension('Time', len(fhrs))
    netcdf_out.createDimension('south_north', len(outlats[:,0]))
    netcdf_out.createDimension('west_east', len(outlons[0,:]))
    netcdf_out.createDimension('sigma', 1)
    times = netcdf_out.createVariable('fhr', int, ('Time'))
    sig = netcdf_out.createVariable('sigma', float, ('sigma'))
    pperfout = netcdf_out.createVariable('practically_perfect', float,
                ('Time', 'south_north', 'west_east'))
    sig[:] = sigma
    times[:] = fhrs

    # Populate outfile with sspf as percentages
    pperfout[:] = sspf_arr*100.
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

def storeNearestNeighborFortran(modelinit, fcsthr, outpath,
                            sixhour=True, obvar="P_HYD",
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
    obvar ----------------- string indicating key of variable in which
                            reliability ob will be stored
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
    copyfile(wrfrefpath, outpath)
    netcdf_out = Dataset(outpath, "a")
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

def storeIdealizedNearestNeighborFortran(ssr_arr, outpath,
                            wrfrefpath='/lustre/research/bancell/aucolema/HWT2016runs'
                                       '/2016050800/wrfoutREFd2',
                            obvar="P_HYD"):
    '''
    Calculate hourly nearest neighbor on WRF grid.
    Overwrite to WRF outfile for fast verification
    with Fortran 77.

    Inputs
    ------
    ssr_arr --------------- 2D array-like of surrogate severe reports
    outpath --------------- filepath of desired earest neighbor output file
                            (WRF reference file will be copied to this path)
    wrfrefpath ------------ path to reference WRF file that has lat/lon info
    obvar ----------------- string indicating key of variable in which
                            ssr_arr will be stored

    Outputs
    -------
    returns NULL, but saves to netCDF outpath.
    '''
    # Create outfile
    copyfile(wrfrefpath, outpath)
    netcdf_out = Dataset(outpath, "a")

    # Populate outfile with gridded observations
    netcdf_out.variables[obvar][0,0,:,:] = ssr_arr[:,:]
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
    def checkSuccess():
        """ Checks reliability input file to make sure arguments were
            ordered correctly.
        """
        try:
            relin = np.genfromtxt("{}reliability.in".format(basedir), dtype=str)
            # Accurate argument order
            args = [str("\'"+probpath+"\'"), str("\'"+obpath+"\'"),
                    str("\'"+outfile+"\'"),
                    fcsthr, str("\'"+variable+"\'"),
                    rthresh, sixhour, nbrhd, rbox_bounds[0],
                    rbox_bounds[1], rbox_bounds[2], rbox_bounds[3]]
            success = True # Assume success initially
            # Ensure that each argument was placed into the proper line of the
            #  reliability input file
            for ind, line in enumerate(relin):
                # If an argument doesn't line up with the rel in arg, set False
                print(str(args[ind]).replace('\\', ''), line)
                if (str(args[ind]).replace('\\', '') != line):
                    success = False
                print(success)
        except:
            success = False
        return success

    def redo_rel_input_file():
        """
            Sometimes function sends arguments into
            reliability input file out of order.
            If this happens, redo.
        """
        if os.path.exists(basedir+'reliability.in'):
            subprocess_cmd("rm {}/reliability.in".format(basedir))
        args = [str("\'"+probpath+"\'"), str("\'"+obpath+"\'"),
                str("\'"+outfile+"\'"),
                fcsthr, str("\'"+variable+"\'"),
                rthresh, sixhour, nbrhd, rbox_bounds[0],
                rbox_bounds[1], rbox_bounds[2], rbox_bounds[3]]
        print("Throwing these into reliability.in")
        print(args)
        with open("{}reliability.in".format(basedir), 'w') as file:
            for arg in args:
                file.write(f"{arg}\n")
        return

    # Get sensitivity input file to define reliability arguments
    esensin = np.genfromtxt(rboxpath)
    rbox_bounds = esensin[4:8] # Response box corners
    # Remove any pre-existing reliability files
    if os.path.exists(basedir+'reliability.in'):
        subprocess_cmd("rm {}/reliability.in".format(basedir))
    # Define arguments
    args = [str("\'"+probpath+"\'"), str("\'"+obpath+"\'"),
            str("\'"+outfile+"\'"),
            fcsthr, str("\'"+variable+"\'"),
            rthresh, sixhour, nbrhd, rbox_bounds[0],
            rbox_bounds[1], rbox_bounds[2], rbox_bounds[3]]
    print("Throwing these into reliability.in")
    print(args)
    # Write to reliability input file
    with open("{}reliability.in".format(basedir), 'w') as file:
        for arg in args:
            file.write(f"{arg}\n")

    # If arguments were mixed up in reliability input file, redo
    # count = 1
    # while checkSuccess() == False:
    #     print("redoing input file, attempt: {}".format(count))
    #     redo_rel_input_file()
    #     try:
    #         relin = np.genfromtxt("{}reliability.in".format(basedir), dtype=str)
    #     except:
    #         redo_rel_input_file()
    #         relin = np.genfromtxt("{}reliability.in".format(basedir), dtype=str)
    #     output_path = relin[2]
    #     print(output_path, "'{}'".format(outfile))
    #     count += 1

    # Remove any pre-existing stored reliability calculations
    if os.path.exists(outfile):
        subprocess_cmd("rm {}".format(outfile))
        print("Removed old reliability output file...")
    print("Running reliability arguments...")
    print("Directory:", package_dir)

    # Yeet that fortran!
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
        endtime_minus_nhrs = len(dat.values)-1 - i
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


def destagger(var, stagger_dim):
    """Return the variable on the unstaggered grid.
    This function destaggers the variable by taking the average of the
    values located on either side of the grid box. Copied from wrf-python
    (https://github.com/NCAR/wrf-python/blob/master/src/wrf/destag.py)

    :param var: `xarray.DataArray`: A variable
            on a staggered grid.
    :param stagger_dim: (:obj:`int`): The dimension index to destagger.
            Negative values can be used to choose dimensions referenced
            from the right hand side (-1 is the rightmost dimension).
    Returns:
        `xarray.DataArray`: The destaggered variable.
    """
    var_shape = var.shape
    num_dims = var.ndim
    stagger_dim_size = var_shape[stagger_dim]

    # Dynamically building the range slices to create the appropriate
    # number of ':'s in the array accessor lists.
    # For example, for a 3D array, the calculation would be
    # result = .5 * (var[:,:,0:stagger_dim_size-2]
    #                    + var[:,:,1:stagger_dim_size-1])
    # for stagger_dim=2.  So, full slices would be used for dims 0 and 1, but
    # dim 2 needs the special slice.
    full_slice = slice(None)
    slice1 = slice(0, stagger_dim_size - 1, 1)
    slice2 = slice(1, stagger_dim_size, 1)

    # default to full slices
    dim_ranges_1 = [full_slice] * num_dims
    dim_ranges_2 = [full_slice] * num_dims

    # for the stagger dim, insert the appropriate slice range
    dim_ranges_1[stagger_dim] = slice1
    dim_ranges_2[stagger_dim] = slice2

    result = .5*(var[tuple(dim_ranges_1)] + var[tuple(dim_ranges_2)])

    return result


def open_wrf_dataset(inname, nest='static', dask=True, chunks=None):
    """
    Runs the WRF Post Processor
    :param inname: string of input file path
    :param nest: string, 'moving' for a moving nest simulation
    :param dask: bool, Specify whether to use dask as backend
    :param chunks: optional dictionary to specify chunks by each dimension. See dask chunks.
    :return: Xarray Dataset of CF-compliant WRF output
    """
    # open the input file
    # specify chunks if not given
    if dask:
        if chunks is None:
            chunks = {'Time': 50, 'west_east': 100, 'west_east_stag': 100,
                      'south_north': 100, 'south_north_stag': 100,
                      'bottom_top': 10, 'bottom_top_stag': 10}

        indata = xr.open_dataset(inname, chunks=chunks)
    else:
        indata = xr.open_dataset(inname, chunks=chunks)

    # Set up output dataset with desired dimensions
    coords = {}
    if nest == 'moving':
        coords['latitude'] = (('time', 'y', 'x'), indata.XLAT.data)
        coords['longitude'] = (('time', 'y', 'x'), indata.XLONG.data)
    else:
        coords['latitude'] = (('y', 'x'), indata.XLAT[0].data)
        coords['longitude'] = (('y', 'x'), indata.XLONG[0].data)
    coords['z'] = (('z'), indata.bottom_top.data)
    coords['time'] = indata.XTIME.data
    ds = xr.Dataset(coords=coords)

    # copy original global attributes and make CF-compliant
    for attr in indata.attrs:
        if attr == 'START_DATE':
            ds.attrs['nest_start_date'] = datetime.strptime(indata.START_DATE,
                                                            '%Y-%m-%d_%H:%M:%S').strftime('%Y-%m-%dT%H:%M:%S')
        elif attr == 'SIMULATION_START_DATE':
            ds.attrs['simulation_start_date'] = datetime.strptime(indata.SIMULATION_START_DATE,
                                                                  '%Y-%m-%d_%H:%M:%S').strftime('%Y-%m-%dT%H:%M:%S')
        elif attr == 'MOAD_CEN_LAT':
            ds.attrs['central_latitude'] = indata.MOAD_CEN_LAT
        # elif attr == 'CEN_LON':
        #     ds.attrs['central_longitude'] = indata.CEN_LON
        elif attr == 'TRUELAT1':
            ds.attrs['true_latitude_1'] = indata.TRUELAT1
        elif attr == 'TRUELAT2':
            ds.attrs['true_latitude_2'] = indata.TRUELAT2
        elif attr == 'STAND_LON':
            ds.attrs['standard_longitude'] = indata.STAND_LON
        elif attr == 'MAP_PROJ_CHAR':
            ds.attrs['projection'] = indata.MAP_PROJ_CHAR
        elif attr == 'WEST-EAST_GRID_DIMENSION':
            continue
        elif attr =='SOUTH-NORTH_GRID_DIMENSION':
            continue
        elif attr == 'BOTTOM-TOP_GRID_DIMENSION':
            continue
        elif attr == 'DX':
            ds.attrs['dx'] = indata.DX
        elif attr == 'DY':
            ds.attrs['dy'] = indata.DY
        else:
            ds.attrs[attr] = indata.attrs[attr]

    # Calculate the model projection x and y coordinates
    # TODO: Add support for more projections
    r = 6370000
    x_model, y_model = lcc_projection(indata, r=r)
    ds['x'] = x_model
    ds.x.attrs['description'] = 'Projection x coordinate'
    ds['y'] = y_model
    ds.y.attrs['description'] = 'Projection y coordinate'

    # add projection information
    ds.attrs['projection'] = 'Lambert Conformal Conic'
    ds.attrs['semimajor_axis'] = r
    ds.attrs['semiminor_axis'] = r
    ds.attrs['ellipse'] = 'sphere'

    # combine and write precipitation variables
    ds['total_precipitation'] = (('time', 'y', 'x'), (indata.RAINC + indata.RAINNC).data)
    ds.total_precipitation.attrs['description'] = 'Total Accumulated Precipitation'
    if indata.RAINC.units == 'mm':
        ds.total_precipitation.attrs['units'] = 'milimeter'
    else:
        warnings.warn('Unknown precipitation unit {} encountered'.format(indata.RAINC.units))
        ds.total_precipitation.attrs['units'] = indata.RAINC.units

    # get surface variables
    # 2-m temperature
    ds['temperature_2m'] = (('time', 'y', 'x'), indata.T2.data)
    ds.temperature_2m.attrs['description'] = 'Temperature at 2 meters'
    if indata.T2.units == 'K':
        ds.temperature_2m.attrs['units'] = 'kelvin'
    else:
        warnings.warn('Unknown temperature unit {} encountered'.format(indata.T2.units))
        ds.temperature_2m.attrs['units'] = indata.T2.units

    # 2-m potential temperature
    ds['potential_temperature_2m'] = (('time', 'y', 'x'), indata.TH2.data)
    ds.potential_temperature_2m.attrs['description'] = 'Potential temperature at 2 meters'
    if indata.TH2.units == 'K':
        ds.potential_temperature_2m.attrs['units'] = 'kelvin'
    else:
        warnings.warn('Unknown temperature unit {} encountered'.format(indata.TH2.units))
        ds.potential_temperature_2m.attrs['units'] = indata.TH2.units

    # 2-m mixing ratio
    ds['mixing_ratio_2m'] = (('time', 'y', 'x'), indata.Q2.data)
    ds.mixing_ratio_2m.attrs['description'] = 'Mixing ratio at 2 meters'
    if indata.Q2.units == 'kg kg-1':
        ds.mixing_ratio_2m.attrs['units'] = 'dimensionless'
    else:
        warnings.warn('Unknown moisture unit {} encountered'.format(indata.Q2.units))
        ds.mixing_ratio_2m.attrs['units'] = indata.Q2.units

    # 10-m winds on model grid
    u10 = indata.U10.data
    v10 = indata.V10.data
    ds['u_10m'] = (('time', 'y', 'x'), u10)
    ds.u_10m.attrs['description'] = 'U-component of wind at 10 meters (model-relative)'
    if indata.U10.units == 'm s-1':
        ds.u_10m.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity units {} encountered'.format(indata.U10.units))
        ds.u_10m.attrs['units'] = indata.U10.units

    ds['v_10m'] = (('time', 'y', 'x'), v10)
    ds.v_10m.attrs['description'] = 'V-component of wind at 10 meters (model-relative)'
    if indata.V10.units == 'm s-1':
        ds.v_10m.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity units {} encountered'.format(indata.V10.units))
        ds.v_10m.attrs['units'] = indata.V10.units

    # 10-m winds earth relative
    sinalpha = indata.SINALPHA.data
    cosalpha = indata.COSALPHA.data
    u10_rot, v10_rot = earth_relative_winds(u10, v10, sinalpha, cosalpha)
    ds['u_10m_earth_relative'] = (('time', 'y', 'x'), u10_rot)
    ds.u_10m_earth_relative.attrs['description'] = 'U-component of wind at 10 meters (earth-relative)'
    if indata.U10.units == 'm s-1':
        ds.u_10m_earth_relative.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity units {} encountered'.format(indata.U10.units))
        ds.u_10m_earth_relative.attrs['units'] = indata.U10.units

    ds['v_10m_earth_relative'] = (('time', 'y', 'x'), v10_rot)
    ds.v_10m_earth_relative.attrs['description'] = 'V-component of wind at 10 meters (earth-relative)'
    if indata.V10.units == 'm s-1':
        ds.v_10m_earth_relative.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity units {} encountered'.format(indata.V10.units))
        ds.v_10m_earth_relative.attrs['units'] = indata.V10.units

    # Surface pressure
    ds['surface_pressure'] = (('time', 'y', 'x'), indata.PSFC.data)
    ds.surface_pressure.attrs['description'] = 'Pressure at surface (not reduced)'
    if indata.PSFC.units == 'Pa':
        ds.surface_pressure.attrs['units'] = 'Pascal'
    else:
        warnings.warn('Unknown pressure unit {} encountered'.format(indata.PSFC.units))
        ds.surface_pressure.attrs['units'] = indata.PSFC.units

    # Get full 3-D variables
    # Height - from geopotential
    hgt = (indata.PH.data + indata.PHB.data) / 9.81
    ds['height_mean_sea_level'] = (('time', 'z', 'y', 'x'), destagger(hgt, 1))
    ds.height_mean_sea_level.attrs['Description'] = 'Height above mean sea level'
    ds.height_mean_sea_level.attrs['units'] = 'meter'

    # Terrain height
    ds['terrain_height'] = (('time', 'y', 'x'), indata.HGT.data)
    ds.terrain_height.attrs['Description'] = 'Height of model terrain'
    if indata.HGT.units == 'm':
        ds.terrain_height.attrs['units'] = 'meter'
    else:
        warnings.warn('Unknown height unit {} encountered'.format(indata.HGT.units))
        ds.terrain_height.attrs['units'] = indata.HGT.units

    # Pressure
    # TODO: find a way to not assume units of input data
    p = indata.P + indata.PB
    ds['pressure'] = (('time', 'z', 'y', 'x'), p.data)
    ds.pressure.attrs['description'] = 'Full model pressure'
    if indata.P.units == 'Pa':
        ds.pressure.attrs['units'] = 'Pascal'
    else:
        warnings.warn('Unknown pressure unit {} encountered'.format(indata.P.units))
        ds.pressure.attrs['units'] = indata.P.units

    theta = indata.T + 300.
    ds['potential_temperature'] = (('time', 'z', 'y', 'x'), theta.data)
    if indata.T.units == 'K':
        ds.potential_temperature.attrs['units'] = 'Kelvin'
    else:
        warnings.warn('Unknown pressure unit {} encountered'.format(indata.T.units))
        ds.potential_temperature.attrs['units'] = indata.T.units

    temp = temperature_from_potential_temperature(ds.pressure, ds.potential_temperature, )
    ds['temperature'] = temp
    if indata.T.units == 'K':
        ds.temperature.attrs['units'] = 'Kelvin'
    else:
        warnings.warn('Unknown pressure unit {} encountered'.format(indata.T.units))
        ds.temperature.attrs['units'] = indata.T.units

    # Vapor mixing ratio
    ds['vapor_mixing_ratio'] = (('time', 'z', 'y', 'x'), indata.QVAPOR.data)
    ds.vapor_mixing_ratio.attrs['description'] = 'Water vapor mixing ratio'
    if indata.QVAPOR.units == 'kg kg-1':
        ds.vapor_mixing_ratio.attrs['units'] = 'dimensionless'
    else:
        warnings.warn('Unknown unit {} encountered'.format(indata.QVAPOR.units))
        ds.vapor_mixing_ratio.attrs['units'] = indata.QVAPOR.units

    # Cloud mixing ratio
    ds['cloud_mixing_ratio'] = (('time', 'z', 'y', 'x'), indata.QCLOUD.data)
    ds.cloud_mixing_ratio.attrs['description'] = 'Cloud water mixing ratio'
    if indata.QCLOUD.units == 'kg kg-1':
        ds.cloud_mixing_ratio.attrs['units'] = 'dimensionless'
    else:
        warnings.warn('Unknown units {} encountered'.format(indata.QCLOUD.units))
        ds.cloud_mixing_ratio.attrs['units'] = indata.QCLOUD.units

    # Rain mixing ratio
    ds['rain_mixing_ratio'] = (('time', 'z', 'y', 'x'), indata.QRAIN.data)
    ds.rain_mixing_ratio.attrs['description'] = 'Rain water mixing ratio'
    if indata.QRAIN.units == 'kg kg-1':
        ds.rain_mixing_ratio.attrs['units'] = 'dimensionless'
    else:
        warnings.warn('Unknown unit {} encountered'.format(indata.QRAIN.units))
        ds.rain_mixing_ratio.attrs['units'] = indata.QRAIN.units

    # Ice mixing ratio
    try:
        ds['ice_mixing_ratio'] = (('time', 'z', 'y', 'x'), indata.QICE.data)
        ds.ice_mixing_ratio.attrs['description'] = 'Ice mixing ratio'
        if indata.QICE.units == 'kg kg-1':
            ds.ice_mixing_ratio.attrs['units'] = 'dimensionless'
        else:
            warnings.warn('Unknown unit {} encountered'.format(indata.QICE.units))
            ds.ice_mixing_ratio.attrs['units'] = indata.QICE.units
    except AttributeError:
        pass

    # Snow mixing ratio
    try:
        ds['snow_mixing_ratio'] = (('time', 'z', 'y', 'x'), indata.QSNOW.data)
        ds.snow_mixing_ratio.attrs['description'] = 'Snow mixing ratio'
        qsnow = indata.QSNOW.data
        if indata.QSNOW.units == 'kg kg-1':
            ds.snow_mixing_ratio.attrs['units'] = 'dimensionless'
        else:
            warnings.warn('Unknown unit {} encountered'.format(indata.QSNOW.units))
            ds.snow_mixing_ration.attrs['units'] = indata.QSNOW.units
    except AttributeError:
        qsnow = None

    # Graupel mixing ratio
    try:
        ds['graupel_mixing_ratio'] = (('time', 'z', 'y', 'x'), indata.QGRAUP.data)
        ds.graupel_mixing_ratio.attrs['description'] = 'Graupel mixing ratio'
        qgraupel = indata.QGRAUP.data
        if indata.QGRAUP.units == 'kg kg-1':
            ds.graupel_mixing_ratio.attrs['units'] = 'dimensionless'
        else:
            warnings.warn('Unknown unit {} encountered'.format(indata.QGRAUP.units))
            ds.graupel_mixing_ratio.attrs['units'] = indata.QGRAUP.units
    except AttributeError:
        qgraupel = None

    # Ice number concentration
    try:
        ds['ice_number_concentration'] = (('time', 'z', 'y', 'x'), indata.QNICE.data)
        ds.ice_number_concentration.attrs['description'] = 'Ice number concentration'
        if indata.QNICE.units == '  kg-1':
            ds.ice_number_concentration.attrs['units'] = 'kiligram**-1'
        else:
            warnings.warn('Unknown unit {} encountered'.format(indata.QNICE.units))
            ds.ice_number_concentration.attrs['units'] = indata.QNICE.units
    except AttributeError:
        pass

    # Rain number concentration
    try:
        ds['rain_number_concentration'] = (('time', 'z', 'y', 'x'), indata.QNRAIN.data)
        ds.rain_number_concentration.attrs['description'] = 'Rain number concentration'
        if indata.QNRAIN.units == '  kg(-1)':
            ds.rain_number_concentration.attrs['units'] = 'kiligram**-1'
        else:
            warnings.warn('Unknown unit {} encountered'.format(indata.QNRAIN.units))
            ds.rain_number_concentration.attrs['units'] = indata.QNRAIN.units
    except AttributeError:
        pass

    # Unstagger the staggered-grid variables
    # u-component of wind model relative
    u = destagger(indata.U.data, -1)
    v = destagger(indata.V.data, -2)

    ds['u_wind'] = (('time', 'z', 'y', 'x'), u)
    ds.u_wind.attrs['description'] = 'U-component of wind (model-relative)'
    if indata.U.units == 'm s-1':
        ds.u_wind.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity unit {} encountered'.format(indata.U.units))
        ds.u_wind.attrs['units'] = indata.U.units

    # v-component of wind
    ds['v_wind'] = (('time', 'z', 'y', 'x'), v)
    ds.v_wind.attrs['description'] = 'V-component of wind (model-relative)'
    if indata.V.units == 'm s-1':
        ds.v_wind.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity unit {} encountered'.format(indata.V.units))
        ds.v_wind.attrs['units'] = indata.V.units

    # earth-relative
    u_rot, v_rot = earth_relative_winds(u, v, sinalpha[:, np.newaxis, ], cosalpha[:, np.newaxis, ])

    ds['u_wind_earth_relative'] = (('time', 'z', 'y', 'x'), u_rot)
    ds.u_wind_earth_relative.attrs['description'] = 'U-component of wind (earth-relative)'
    if indata.U.units == 'm s-1':
        ds.u_wind_earth_relative.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity unit {} encountered'.format(indata.U.units))
        ds.u_wind_earth_relative.attrs['units'] = indata.U.units

    # v-component of wind
    ds['v_wind_earth_relative'] = (('time', 'z', 'y', 'x'), v_rot)
    ds.v_wind_earth_relative.attrs['description'] = 'V-component of wind (earth-relative)'
    if indata.V.units == 'm s-1':
        ds.v_wind_earth_relative.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity unit {} encountered'.format(indata.V.units))
        ds.v_wind_earth_relative.attrs['units'] = indata.V.units

    # w-component of wind
    ds['w_wind'] = (('time', 'z', 'y', 'x'), destagger(indata.W.data, -3))
    ds.w_wind.attrs['description'] = 'W-component of wind'
    if indata.W.units == 'm s-1':
        ds.w_wind.attrs['units'] = 'meter/second'
    else:
        warnings.warn('Unknown velocity unit {} encountered'.format(indata.W.units))
        ds.w_wind.attrs['units'] = indata.W.units

    # COSALPHA and SINALPHA for model grid to earth-relative rotation
    ds['cosalpha'] = (('time', 'y', 'x'), sinalpha)
    ds.cosalpha.attrs['description'] = 'Cosine Alpha term for earth-relative grid rotation'

    ds['sinalpha'] = (('time', 'y', 'x'), cosalpha)
    ds.sinalpha.attrs['description'] = 'Sine Alpha term for earth-relative grid rotation'

    # TKE from PBL scheme
    try:
        ds['tke'] = (('time', 'z', 'y', 'x'), destagger(indata.TKE_PBL.data, -3))
        ds.tke.attrs['description'] = 'Turbulence Kinetic Energy (TKE) from PBL scheme'
        if indata.TKE_PBL.units == 'm2 s-2':
            ds.w_wind.attrs['units'] = 'meter**2/second**2'
        else:
            warnings.warn('Unknown unit {} encountered'.format(indata.TKE_PBL.units))
            ds.tke.attrs['units'] = indata.TKE_PBL.units
    except AttributeError:
        pass

    refl = simulated_reflectivity(ds['pressure'].data, ds.temperature.data, ds['vapor_mixing_ratio'].data,
                                  ds['rain_mixing_ratio'].data, snow_mixing_ratio=qsnow,
                                  graupel_mixing_ratio=qgraupel)
    ds['simulated_reflectivity'] = (('time', 'z', 'y', 'x'), refl)
    ds.simulated_reflectivity.attrs['units'] = 'dBZ'
    return ds


def get_nearest(src_points, candidates, k_neighbors=1):
    """Find nearest neighbors for all source points from a set of candidate points"""

    # Create tree from the candidate points
    tree = BallTree(candidates, leaf_size=15, metric='haversine')
    distances, indices = tree.query(src_points, k=k_neighbors)

    # Transpose to get distances and indices into arrays
    distances = distances.transpose()
    indices = indices.transpose()

    # Get closest indices and distances (i.e. array at index 0)
    # note: for the second closest points, you would take index 1, etc.
    closest = indices[0]
    closest_dist = distances[0]

    # Return indices and distances
    return closest, closest_dist


def nearest_lat_lon_index(point, latitudes, longitudes):
    """
    Find the index values of the nearest grid point to a desired point.
    The Haversine distance formula is used, and inputs of lat/lon pairs are
    required.

    :param point: list, list of latitude, longitude for desired nearest point
    :param latitudes: dask array of 2-D latitude grid
    :param longitudes: dask array of 2-D longitude grid
    :return: index values of nearest point on grid
    """
    stacked = da.stack((latitudes.flatten(), longitudes.flatten())).transpose()
    closest = get_nearest(da.array(point).reshape(-1, 1).transpose(), stacked)
    idx = np.unravel_index(closest[0], latitudes.shape)
    return idx


def earth_relative_winds(u, v, sinalpha, cosalpha):
    """
    Rotate model-relative wind components to earth-relative

    :param u: x-wind component (model relative)
    :param v: y-wind component (model relative)
    :param sinalpha:
    :param cosalpha:
    :return: u_rot, v_rot: u and v wind components rotated to earth relative
    """
    u_rot = u * cosalpha - v * sinalpha
    v_rot = v * cosalpha + u * sinalpha
    return u_rot, v_rot


def lcc_projection(indata, r=6370000):
    """
    Define projection coordinates for WRF Lambert Conformal Conic grid
    :param indata: Xarray.Dataset containing WRF output as netCDF
    :param r: radius of earth (spherical)
    :return: x, y: x and y coordinate arrays
    """
    wrf_proj = pyproj.Proj(proj='lcc',  # projection type: Lambert Conformal Conic
                           lat_1=indata.TRUELAT1, lat_2=indata.TRUELAT2,  # Cone intersects with the sphere
                           lat_0=indata.MOAD_CEN_LAT, lon_0=indata.STAND_LON,  # Center point
                           a=r, b=r)  # This is it! The Earth is a perfect sphere
    wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
    e, n = pyproj.transform(wgs_proj, wrf_proj, indata.CEN_LON, indata.CEN_LAT)
    # Grid parameters
    dx, dy = indata.DX, indata.DY
    nx, ny = float(indata.dims['west_east']), float(indata.dims['south_north'])
    # Down left corner of the domain
    x0 = -(nx - 1) / 2. * dx + e
    y0 = -(ny - 1) / 2. * dy + n

    x = da.arange(nx) * dx + x0
    y = da.arange(ny) * dy + y0
    return x, y


def simulated_reflectivity(pressure, temperature, vapor_mixing_ratio, liquid_mixing_ratio, snow_mixing_ratio=None,
                           graupel_mixing_ratio=None, use_varint=False, use_liqskin=False):
    """
    Calculate the simulated reflectivity factor from model output.
    Ported from RIP fortran calculation used in the WRF-Python package.
    :param pressure: model pressure in Pa
    :param temperature: model temperature in Kelvin
    :param vapor_mixing_ratio: water vapor mixing ratio
    :param liquid_mixing_ratio: liquid water mixing ratio
    :param snow_mixing_ratio: snow mixing ratio, optional
    :param graupel_mixing_ratio: graupel mixing ratio, optional
    :param use_varint: When set to False,
        the intercept parameters are assumed constant
        (as in MM5's Reisner-2 bulk microphysical scheme).
        When set to True, the variable intercept
        parameters are used as in the more recent version of Reisner-2
        (based on Thompson, Rasmussen, and Manning, 2004, Monthly weather
        Review, Vol. 132, No. 2, pp. 519-542.).
    :param use_liqskin: When set to True, frozen particles
        that are at a temperature above freezing are assumed to scatter
        as a liquid particle.  Set to False to disable.

    This routine computes equivalent reflectivity factor (in dBZ) at
    each model grid point.  In calculating Ze, the RIP algorithm makes
    assumptions consistent with those made in an early version
    (ca. 1996) of the bulk mixed-phase microphysical scheme in the MM5
    model (i.e., the scheme known as "Resiner-2").  For each species:

    1. Particles are assumed to be spheres of constant density.  The
    densities of rain drops, snow particles, and graupel particles are
    taken to be rho_r = rho_l = 1000 kg m^-3, rho_s = 100 kg m^-3, and
    rho_g = 400 kg m^-3, respectively. (l refers to the density of
    liquid water.)

    2. The size distribution (in terms of the actual diameter of the
    particles, rather than the melted diameter or the equivalent solid
    ice sphere diameter) is assumed to follow an exponential
    distribution of the form N(D) = N_0 * exp( lambda*D ).

    3. If ivarint=0, the intercept parameters are assumed constant
    (as in early Reisner-2), with values of 8x10^6, 2x10^7,
    and 4x10^6 m^-4, for rain, snow, and graupel, respectively.
    If ivarint=1, variable intercept parameters are used, as
    calculated in Thompson, Rasmussen, and Manning (2004, Monthly
    Weather Review, Vol. 132, No. 2, pp. 519-542.)

    4. If iliqskin=1, frozen particles that are at a temperature above
    freezing are assumed to scatter as a liquid particle.

    More information on the derivation of simulated reflectivity in
    RIP can be found in Stoelinga (2005, unpublished write-up).
    Contact Mark Stoelinga (stoeling@atmos.washington.edu) for a copy.
    """
    # Set values for constants with variable intercept
    R1 = 1e-15
    RON = 8e6
    RON2 = 1e10
    SON = 2e7
    GON = 5e7
    RON_MIN = 8e6
    RON_QR0 = 0.00010
    RON_DELQR0 = 0.25*RON_QR0
    RON_CONST1R = (RON2-RON_MIN)*0.5
    RON_CONST2R = (RON2+RON_MIN)*0.5

    # set constant intercepts
    rno_l = 8e6
    rno_s = 2e7
    rno_g = 4e6

    qvapor = da.clip(vapor_mixing_ratio, 0., None)
    qliquid = da.clip(liquid_mixing_ratio, 0., None)

    # If qgraupel but not qsnow, set qgraupel = qsnow
    if snow_mixing_ratio is None:
        if graupel_mixing_ratio is None:
            qsnow = da.zeros_like(qliquid)
            qgraupel = da.zeros_like(qliquid)
        else:
            qgraupel = da.clip(graupel_mixing_ratio, 0., None)
            qsnow = da.zeros_like(graupel_mixing_ratio)
            qsnow[temperature <= 273.15] = qgraupel[temperature <= 273.15]
    else:
        qsnow = da.clip(snow_mixing_ratio, 0., None)
        qgraupel = da.clip(graupel_mixing_ratio, 0., None)

    # density for liquid, snow, and graupel (kg m-3)
    rho_l = 1000.  # liquid
    rho_i = 100.  # snow
    rho_g = 400.  # graupel

    # constant evaluation of gamma distribution
    gamma = 720.

    # Alpha constant
    alpha = 0.224

    # constant multiplication factors
    factor_l = gamma * 1e18 * (1./(np.pi*rho_l))**1.75
    s = gamma * 1e18 * (1./(np.pi*rho_i))**1.75 * (rho_i/rho_l)**2 * alpha
    g = gamma * 1e18 * (1./(np.pi*rho_g))**1.75 * (rho_g/rho_l)**2 * alpha

    # calculate virtual temperature
    virtual_t = virtual_temperature(temperature, qvapor)

    # dry gas constant
    Rd = 287.
    rho_air = pressure/(Rd*virtual_t)

    # adjust for brightband if use_liqskin=True
    if use_liqskin:
        raise NotImplementedError('Liquid skin correction not implemented')
        # factor_s = da.full_like(temperature, s)
        # factor_g = da.full_like(temperature, g)
        # try:
        #     factor_s[temperature >= 273.15] = factor_s[temperature >= 273.15] / da.array([alpha])
        #     factor_g[temperature >= 273.15] = factor_g[temperature >= 273.15] / da.array([alpha])
        # except ValueError:
        #     factor_s = s
        #     factor_g = g
    else:
        factor_s = s
        factor_g = g

    # calculate variable intercept if use_varint=True
    if use_varint:
        raise NotImplementedError('Variable intercepts not yet implemented')
        # temp_c = da.clip(temperature-273.15, temperature.min(), -0.001)
        # sonv = MIN(2.0D8, 2.0D6*EXP(-0.12D0*temp_c))
        #
        # gonv = gon
        # IF (qgr(i,j,k) .GT. R1) THEN
        #     gonv = 2.38D0 * (PI*RHO_G/(rhoair*qgr(i,j,k)))**0.92D0
        #     gonv = MAX(1.D4, MIN(gonv,GON))
        # END IF
        #
        # ronv = RON2
        # IF (qra(i,j,k) .GT. R1) THEN
        #     ronv = RON_CONST1R*TANH((RON_QR0 - qra(i,j,k))/RON_DELQR0) + RON_CONST2R
        # END IF
    else:
        ronv = rno_l
        sonv = rno_s
        gonv = rno_g

    # Total equivalent reflectivity factor (z_e, in mm^6 m^-3) is
    # the sum of z_e for each hydrometeor species:
    z_e = (((factor_l*(rho_air*qliquid)**1.75)/(ronv**.75)) +
           ((factor_s*(rho_air*qsnow)**1.75)/(sonv**.75)) +
           ((factor_g*(rho_air*qgraupel)**1.75)/(gonv**.75)))

    # Adjust small values of Z_e so that dBZ is no lower than -30
    z_e = da.clip(z_e, .001, None)

    # Convert to dBZ
    dbz = 10.*da.log10(z_e)
    return dbz


def virtual_temperature(temperature, mixing, molecular_weight_ratio=0.622):
    r"""Calculate virtual temperature.

    This calculation must be given an air parcel's temperature and mixing ratio.
    The implementation uses the formula outlined in [Hobbs2006]_ pg.80. Taken from metpy.calc
    and modified for Dask support.

    Parameters
    ----------
    temperature:
        air temperature
    mixing :
        dimensionless mass mixing ratio
    molecular_weight_ratio : float, optional
        The ratio of the molecular weight of the constituent gas to that assumed
        for air. Defaults to the ratio for water vapor to dry air.
        (:math:`\epsilon\approx0.622`).

    Returns
    -------
        The corresponding virtual temperature of the parcel

    Notes
    -----
    .. math:: T_v = T \frac{\text{w} + \epsilon}{\epsilon\,(1 + \text{w})}

    """
    return temperature * ((mixing + molecular_weight_ratio)
                          / (molecular_weight_ratio * (1 + mixing)))


def exner_function(pressure, reference_pressure=P0):
    r"""Calculate the Exner function. From metpy.calc.
    .. math:: \Pi = \left( \frac{p}{p_0} \right)^\kappa
    This can be used to calculate potential temperature from temperature (and visa-versa),
    since
    .. math:: \Pi = \frac{T}{\theta}
    Parameters
    ----------
    pressure :
        total atmospheric pressure, units should match reference pressure (default is Pa)
    reference_pressure :
        The reference pressure against which to calculate the Exner function, defaults to
        metpy.constants.P0
    Returns
    -------
        The value of the Exner function at the given pressure
    See Also
    --------
    potential_temperature
    temperature_from_potential_temperature
    """
    return (pressure / reference_pressure)**kappa


def temperature_from_potential_temperature(pressure, potential_temperature, reference_pressure=P0):
    r"""Calculate the temperature from a given potential temperature.
    Uses the inverse of the Poisson equation to calculate the temperature from a
    given potential temperature at a specific pressure level. Taken from metpy.calc and modified for Dask.
    Parameters
    ----------
    pressure :
        total atmospheric pressure, units should match reference pressure (default is Pa).
    potential_temperature :
        potential temperature
    reference_pressure :
        The reference pressure against which to calculate the Exner function, defaults to
        metpy.constants.P0
    Returns
    -------
        The temperature corresponding to the potential temperature and pressure.
    See Also
    --------
    dry_lapse
    potential_temperature
    Notes
    -----
    Formula:
    .. math:: T = \Theta (P / P_0)^\kappa
    """
    return potential_temperature * exner_function(pressure, reference_pressure=reference_pressure)
