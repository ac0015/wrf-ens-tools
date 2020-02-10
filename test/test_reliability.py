import numpy as np; import os
from wrf_ens_tools.post import storeReliabilityRboxFortran
from datetime import datetime
from netCDF4 import Dataset

# Create dummy forecast data
xdim, ydim, zdim = 588, 540, 40
fake_fcst_data = np.random.uniform(low=0., high=100., size=(1, zdim, ydim, xdim))
fhr = 28

# Create dummy observational data
prob_hit = 0.3
fake_ob_data = np.random.choice([0, 1], p=[1-prob_hit, prob_hit], size=(1, ydim, xdim))

# Pull real lat/lons
real_probs = "wrfoutREFd2"
latlondat = Dataset(real_probs)
lats = latlondat.variables['XLAT'][:]
lons = latlondat.variables['XLONG'][:]

# Store dummy data
dummy_probpath = "dummy_wrfoutprobs.prob"
dummy_dat = Dataset(dummy_probpath, 'w')
dummy_dat.DX = 4000.
dummy_dat.createDimension('Time', size=None)
dummy_dat.createDimension('bottom_top', size=zdim)
dummy_dat.createDimension('south_north', size=ydim)
dummy_dat.createDimension('west_east', size=xdim)
dummy_dat.createVariable('P_HYD', float, dimensions=('Time', 'bottom_top',
                                                    'south_north', 'west_east'))
dummy_dat['P_HYD'][:] = fake_fcst_data
dummy_dat.createVariable('XLONG', float, dimensions=('Time', 'south_north',
                                                    'west_east'))
dummy_dat['XLONG'][:] = lons
fcsthr = dummy_dat.createVariable('fhr', int, dimensions=('Time'))
fcsthr[:] = fhr
dummy_dat.createVariable('XLAT', float, dimensions=('Time', 'south_north',
                                                    'west_east'))
dummy_dat['XLAT'][:] = lats
dummy_dat.close()

# Copy dummy_dat into ob dat
dummy_obpath = "dummy_pperf.prob"
dummy_dat = Dataset(dummy_obpath, 'w')
dummy_dat.createDimension('Time', size=None)
dummy_dat.createDimension('bottom_top', size=zdim)
dummy_dat.createDimension('south_north', size=ydim)
dummy_dat.createDimension('west_east', size=xdim)
dummy_dat.createVariable('nearest_neighbor', float, dimensions=('Time',
                                                    'south_north', 'west_east'))
dummy_dat['nearest_neighbor'][:] = fake_ob_data
dummy_dat.createVariable('XLONG', float, dimensions=('Time', 'south_north',
                                                    'west_east'))
dummy_dat['XLONG'][:] = lons
dummy_dat.createVariable('XLAT', float, dimensions=('Time', 'south_north',
                                                    'west_east'))
fcsthr = dummy_dat.createVariable('fhr', int, dimensions=('Time'))
fcsthr[:] = fhr
dummy_dat['XLAT'][:] = lats
dummy_dat.close()

# Create fake run for testing
run = datetime(2018, 5, 28, 0)
rboxpath = "esens.in"
sixhr = True

# def test_scipy_reliability_rbox():
#     # Test reliability over response box
#     bins, fcst_freq_rbox, ob_hr_rbox = scipyReliabilityRbox(dummy_probpath, run, fhr,
#                     obpath=dummy_obpath, var='updraft_helicity',
#                     thresh=100., rboxpath=rboxpath, sixhr=sixhr, nbrhd=0.)
#     assert(((ob_hr_rbox - prob_hit) < 0.05).all() == True)
