import numpy as np
from netCDF4 import Dataset
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
import coordinateSystems as cs
from scipy import interpolate
from siphon import ncss
import subprocess

###########################################################
# Austin Coleman
# 12/27/2017
#
# Python library to interpolate external analyses from
# desired date/time to a desired WRF grid and save
# off the interpolated analysis to a netCDF4 file.
##########################################################

# For more convenient subprocess runs (courtesy of Stack Overflow)
def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)
    return

# Clunky alternative to strftime. TO-DO: Move all calls of this to strftime()
def fromDatetime(date, interp=False):
    """
    Slices components of a datetime.datetime object and returns
    as a tuple of strings to be compatible with sensitivity and
    RAP file naming conventions.
    """
    year = str(date.year)
    # Convert any 1-digit integers to 2-digit strings
    if len(str(date.month)) == 1: month = "0"+str(date.month)
    else: month = str(date.month)
    if len(str(date.day)) == 1: day = "0"+str(date.day)
    else: day = str(date.day)
    if len(str(date.hour)) == 1: hour = "0"+str(date.hour)
    else: hour = str(date.hour)
    if interp: hour += "00" # Append two zeroes for 4-digit hr
    return year, month, day, hour

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
    interp = interpolate.LinearNDInterpolator(coords_from, Z, fill_value=9e9)
    # Interpolate to new grid
    interpolated_z = interp(grid2y, grid2x)
    return interpolated_z

def interpRAPtoWRF(yr, mo, day, hr, wrfref):
    """
    Method to interpolate a RAP analysis pulled from the NCSS server
    to a user-provided WRF grid for a specific day and time.

    Inputs
    ------
    yr ----- 4-digit year as string
    month -- 2-digit month (i.e. 02 for February) as string
    day ---- 2-digit day as string
    hr ----- 4-digit hour (i.e. 1200 for 12 UTC) as string
    wrfref - filepath to wrfout grid file to which RAP will be interpolated

    Outputs
    -------
    returns NULL but will produce a netCDF4 file with interpolated RAP in the
    directory where method was called. Also produces a test ".png" file comparing
    RAP to interpolated RAP to ensure that the interpolation worked.
    """
    # Read RAP datset
    rap = ncss.NCSS("https://www.ncei.noaa.gov/thredds/ncss/grid/rap130anl/{}{}/{}{}{}/rap_130_{}{}{}_{}_000.grb2".format(yr,
                    mo, yr, mo, day, yr, mo, day, hr))

    # Develop query
    query = ncss.NCSSQuery()
    query.variables('all')
    query.add_lonlat()

    # Pull RAP data
    rap_anl = rap.get_data(query)
    rlons, rlats = rap_anl.variables['lon'][:,:], rap_anl.variables['lat'][:,:]

    # Start pulling variables
    gph = rap_anl.variables['Geopotential_height_isobaric']
    isobaric = rap_anl.variables['isobaric']
    temp = rap_anl.variables['Temperature_isobaric']
    uwind = rap_anl.variables['u-component_of_wind_isobaric']
    vwind = rap_anl.variables['v-component_of_wind_isobaric']
    hgt_above_gnd = rap_anl.variables['height_above_ground']
    hgt_above_gnd4 = rap_anl.variables['height_above_ground4']
    temp_hgt = rap_anl.variables['Temperature_height_above_ground']
    dwpt_hgt = rap_anl.variables['Dewpoint_temperature_height_above_ground']
    q_hgt = rap_anl.variables['Specific_humidity_height_above_ground']
    slp = rap_anl.variables['MSLP_MAPS_System_Reduction_msl'][0,:,:]
    uwind_hgt = rap_anl.variables['u-component_of_wind_height_above_ground']
    vwind_hgt = rap_anl.variables['v-component_of_wind_height_above_ground']

    # Designate vertical levels of interest
    lev925 = np.where(isobaric[:] == 92500.)[0][0]
    lev850 = np.where(isobaric[:] == 85000.)[0][0]
    lev700 = np.where(isobaric[:] == 70000.)[0][0]
    lev500 = np.where(isobaric[:] == 50000.)[0][0]
    lev300 = np.where(isobaric[:] == 30000.)[0][0]
    lev2m = np.where(hgt_above_gnd[:] == 2.)[0][0]
    lev10m = np.where(hgt_above_gnd4[:] == 10.)[0][0]

    # Specific variables
    gph300 = gph[0,lev300,:,:]
    gph500 = gph[0,lev500,:,:]
    gph700 = gph[0,lev700,:,:]
    gph850 = gph[0,lev850,:,:]
    temp300 = temp[0,lev300,:,:]
    temp500 = temp[0,lev500,:,:]
    temp700 = temp[0,lev700,:,:]
    temp850 = temp[0,lev850,:,:]
    temp925 = temp[0,lev925,:,:]
    u300 = uwind[0,lev300,:,:]
    u500 = uwind[0,lev500,:,:]
    u700 = uwind[0,lev700,:,:]
    u850 = uwind[0,lev850,:,:]
    u925 = uwind[0,lev925,:,:]
    v300 = vwind[0,lev300,:,:]
    v500 = vwind[0,lev500,:,:]
    v700 = vwind[0,lev700,:,:]
    v850 = vwind[0,lev850,:,:]
    v925 = vwind[0,lev925,:,:]
    t2 = temp_hgt[0,lev2m,:,:]
    q2 = q_hgt[0,lev2m,:,:]
    td2 = dwpt_hgt[0,lev2m,:,:]
    u10 = uwind_hgt[0,lev10m,:,:]
    v10 = vwind_hgt[0,lev10m,:,:]

    # Variables to interpolate
    sensvarlist = [gph300,gph500,gph700,gph850,temp300,temp500,temp700,temp850,temp925,
                  u300,u500,u700,u850,u925,v300,v500,v700,v850,v925,slp,t2,q2,td2,u10,v10]
    # TO-DO: Change sensitivity strings to get rid of leading numeric values.
    #        Makes formatting in the netCDF file really gross-looking.
    sensstringslist = ["300 hPa GPH","500 hPa GPH","700 hPa GPH","850 hPa GPH","300 hPa T","500 hPa T",
                       "700 hPa T","850 hPa T","925 hPa T","300 hPa U-Wind","500 hPa U-Wind",
                       "700 hPa U-Wind","850 hPa U-Wind","925 hPa U-Wind","300 hPa V-Wind","500 hPa V-Wind",
                       "700 hPa V-Wind","850 hPa V-Wind","925 hPa V-Wind","SLP","2m Temp","2m Q",
                       "2m Dewpt","10m U-Wind","10m V-Wind"]

    # Initialize geographic coordinate system
    geo = cs.GeographicSystem()

    # Read one of our WRF files in for grid purposes
    wrf_d1 = Dataset(wrfref)
    lons, lats = wrf_d1.variables['XLONG'][0], wrf_d1.variables['XLAT'][0]
    wrf_idim = len(lons[0,:])
    wrf_jdim = len(lats[:,0])

    # Convert to Earth-Center-Earth-Fixed for common origin
    wecefx, wecefy, wecefz = geo.toECEF(lons, lats, lats)
    recefx, recefy, recefz = geo.toECEF(rlons, rlats, rlats)

    # Write interpolated variables to netCDF
    interpnc = Dataset(str("RAP_interp_to_WRF_{}{}{}{}.nc".format(yr,
                           mo, day, hr)), "w", format="NETCDF4")
    interpnc.createDimension('lat', wrf_jdim)
    interpnc.createDimension('lon', wrf_idim)
    interpnc.createDimension('time', None)
    xlat = interpnc.createVariable("XLAT", float, dimensions=('lat','lon'))
    xlat[:,:] = lats
    xlon = interpnc.createVariable("XLONG", float, dimensions=('lat','lon'))
    xlon[:,:] = lons

    # Interpolate and save!!
    RAPinterp = np.zeros((len(sensvarlist), len(lats[:,0]), len(lons[0,:])))
    for i in range(len(sensvarlist)):
        RAPinterp[i] = bilinear_interp(recefx, recefy, wecefx, wecefy, sensvarlist[i])
        var = interpnc.createVariable(sensstringslist[i].replace(" ","_"),
                                RAPinterp[i].dtype, dimensions=('lat','lon'))
        var[:,:] = RAPinterp[i]
    interpnc.close()
    return

def plotRAPinterp(rapfile, var="500_hPa_GPH"):
    """
    Function mainly meant to sanity-check that a RAP
    interpolation worked correctly.

    Inputs
    ------
    rapfile --- absolute filepath of desired
                interpolated RAP file to test.
    var ------- variable string of the variable
                to plot. Defaults to "500500_hPa_GPH".

    Outputs
    -------
    returns NULL but saves an image called
    "RAP_Analysis_{var}.png"
    """
    rap = Dataset(rapfile)
    lons, lats = rap.variables['XLONG'], rap.variables['XLAT']
    rapvar = rap.variables[var]
    rapvar_masked = np.ma.masked_array(rapvar, mask=(rapvar[:,:]==9e9))
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    state_borders = cfeat.NaturalEarthFeature(category='cultural',
               name='admin_1_states_provinces_lakes', scale='50m', facecolor='None')
    ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
    ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
    ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
    ax.set_extent([lons[0,:].min(), lons[0,:].max(), lats[:,0].min(), lats[:,0].max()])
    # Specific 500 GPH levels to compare to SPC Mesoanalysis contours (sanity check)
    if var == "500_hPa_GPH":
       clevels = np.arange(5280,5940,60)
       interp = ax.contour(lons[:,:], lats[:,:], rapvar_masked[:,:], clevels, transform=ccrs.PlateCarree(),
                       alpha=0.7, antialiased=True, colors='k', linewidths=2.)
       ax.clabel(interp)
    else:
       cflevels = np.linspace(np.min(rapvar_masked[:,:]), np.max(rapvar_masked[:,:]), 21)
       interp = ax.contourf(lons[:,:], lats[:,:], rapvar_masked[:,:], cflevels, transform=ccrs.PlateCarree(),
                       alpha=0.7, antialiased=True)
       fig.colorbar(interp, fraction=0.046, pad=0.04, orientation='horizontal')
    ax.set_title(var.replace("_"," ")+" RAP Analysis")
    plt.savefig("RAP_Analysis_"+var)
    plt.close()
    return
