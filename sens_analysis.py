import numpy as np
import sys
from netCDF4 import Dataset
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from scipy.ndimage.filters import gaussian_filter

##########################################################################
# Script to plot ensemble sensitivity for all specified variables
#  via geospatial and scatterplot analyses.
#
# Slice information
# ----------------
# smat/mmat (sensitivity vars and means)
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
#########################################################################

# Pull in sensitivity data and user input
wrfsens = Dataset('wrfout.sens') # wrfout from esensSPC
meta = np.genfromtxt('esens.in') # input for esensSPC
senstime = int(meta[1])
rtime = int(meta[2])
rfunc = int(meta[3])
slon = meta[4]
elon = meta[5]
slat = meta[6]
elat = meta[7]

# Define response function description and units
rfuncstr = ''
runits = ''
if rfunc == 1:
    rfuncstr = 'Avg Sim Refl'
    runits = 'dBZ'
elif rfunc == 2:
    rfuncstr = 'Max Sim Refl'
    runits = 'dBZ'
elif rfunc == 3:
    rfuncstr = 'Avg UH'
    runits = r'$\frac{m^2}{s^2}$'
elif rfunc == 4:
    rfuncstr = 'Max UH'
    runits = r'$\frac{m^2}{s^2}$'
elif rfunc == 5:
    rfuncstr = 'Accum PCP'
    runits = 'Inch'
elif rfunc == 6:
    rfuncstr = 'Avg 10m Wind Speed'
    runits = r'\frac{m}{s}'

# Define sensitivity variable labels
sens_vars = ['300 hPa GPH', '500 hPa GPH', '700 hPa GPH', '850 hPa GPH', '300 hPa Temp (C)', '500 hPa Temp (C)', 
             '700 hPa Temp (C)', '850 hPa Temp (C)', '925 hPa Temp (C)', '300 hPa U-Wind (m per s)',
             '500 hPa U-Wind (m per s)', '700 hPa U-Wind (m per s)', '850 hPa U-Wind (m per s)',
             '925 hPa U-Wind (m per s)', '300 hPa V-Wind (m per s)', '500 hPa V-Wind (m per s)',
             '700 hPa V-Wind (m per s)', '850 hPa V-Wind (m per s)', '925 hPa V-Wind (m per s)',
             '850 hPa Q', 'SLP', '2 m Temp', '2 m Q', '2 m Dewpt', '10 m U-Wind', '10 m V-Wind']

# Pull sensitivty and mean matrices from wrfsens
smat = wrfsens.variables['P_HYD'][0]
mmat = wrfsens.variables['QICE'][0]

# Define response function box
lons, lats = wrfsens.variables['XLONG'][0], wrfsens.variables['XLAT'][0]
width = elon - slon
height = elat - slat

# Mask underground and missing values
mmask = (mmat>=9e9)
mmask2 = (mmat<=-9999.0)
smask = (smat>=9e9)
mmat_masked = np.ma.masked_array(data=mmat, mask=mmask)
mmat_masked2 = np.ma.masked_array(data=mmat_masked, mask=mmask2)
smat_masked = np.ma.masked_array(data=smat, mask=smask)

# Set up max sens arrays
sens_max = np.zeros((len(sens_vars)))

for i in range(len(sens_vars)):
    print(sens_vars[i])
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    rbox = patches.Rectangle((slon, slat), width, height, transform=ccrs.PlateCarree(), fill=False, color='green', linewidth=2., 
                             zorder=3.)
    state_borders = cfeat.NaturalEarthFeature(category='cultural',
        name='admin_1_states_provinces_lakes', scale='50m', facecolor='None') 
    ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
    ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
    ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
    ax.set_extent([lons[0,:].min(), lons[0,:].max(), lats[:,0].min(), lats[:,0].max()])
    step = int((np.max(mmat_masked2[i,:,:]) - np.min(mmat_masked2[i,:,:]))*0.15)
    if step == 0:
        step += 1
    clevels = np.arange(int(np.min(mmat_masked2[i,:,:])), int(np.max(mmat_masked2[i,:,:])), step)
    cflevels = np.linspace(-50., 50., 101)
    sensmean = ax.contour(lons, lats, gaussian_filter(mmat_masked2[i,:,:], 1), clevels, transform=ccrs.PlateCarree(), colors='black')
    esens = ax.contourf(lons, lats, smat_masked[i,:,:], cflevels, transform=ccrs.PlateCarree(), cmap='RdBu_r', 
                             alpha=0.7, extend='both', antialiased=True)
    #p = ax.scatter(lons[187, 123], lats[187, 123], transform=ccrs.PlateCarree(), s=100., c='orange', zorder=12)
    #print(smat_masked[i,187, 123])
    esens.cmap.set_under(color=u'navy')
    esens.cmap.set_over(color=u'darkred')
    ax.add_patch(rbox)
    cbar = fig.colorbar(esens, fraction=0.046, pad=0.04, orientation='horizontal', label='Sensitivity')
    ax.clabel(sensmean, fmt='%4.0f')
    ax.set_title('Sensitivity of '+rfuncstr+' at f'+str(rtime)+' to '+sens_vars[i]+' at f'+str(senstime))
    plt.savefig(rfuncstr.replace(" ", "")+'f'+str(rtime)+'_to_'+sens_vars[i].replace(" ","")+'f'+str(senstime)) 
    plt.close()   
