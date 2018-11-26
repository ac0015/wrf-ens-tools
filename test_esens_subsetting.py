import pytest
import esens_subsetting
import numpy as np
from netCDF4 import Dataset

# Dimensions of test array
xdim, ydim, zdim = 300, 300, 40
ensnum = 10
sensvar, sensind = "500_hPa_GPH", 1
memvalvar, analysisvar = "GPH_500", "500_hPa_GPH"
# Filepaths for storing dummy data
senspath = 'test/sens_rand.nc'
analysispath = 'test/fake_RAP_interp_analysis.nc'
memvalspath = 'test/mem_vals.nc'

# Define arbitrary sensitivity field
sens_field = np.random.rand(1, ydim, xdim)
# "Truth" will consist of an array of ones
analysis = np.ones_like(sens_field)
# Members will all have values of 2
mem_sens_vals = np.ones((ensnum, ydim, xdim)) * 2
# Make member 5 the best member
mem_sens_vals[4,:,:] = 1.25
# Known error field - should all be ones aside from member 5
error_field = np.zeros_like(mem_sens_vals)
sensweighted_error_field = np.zeros_like(error_field)
for k in range(ensnum):
    error_field[k] = mem_sens_vals[k] - analysis
    sensweighted_error_field[k] = sens_field[0] * error_field[k]
print("Sanity check - Does max of error fields equal one?",
    np.max(error_field))
# Store arrays to netCDF4 that emulate a wrfout.sens file
sensdat = Dataset(senspath, 'w')
sensdat.createDimension('Time', size=None)
sensdat.createDimension('bottom_top', size=zdim)
sensdat.createDimension('south_north', size=ydim)
sensdat.createDimension('west_east', size=xdim)
sensdat.createVariable('P_HYD', float, dimensions=('Time', 'bottom_top',
                                                    'south_north', 'west_east'))
sensdat.createVariable('XLONG', float, dimensions=('Time', 'south_north',
                                                    'west_east'))
sensdat['XLONG'][:] = np.ones_like(sens_field)
sensdat.createVariable('XLAT', float, dimensions=('Time', 'south_north',
                                                    'west_east'))
sensdat['XLAT'][:] = np.ones_like(sens_field)

sensdat['P_HYD'][:,sensind,:,:] = sens_field
sensdat.close()
# Store analysis data to netCDF4 using the RAP interpolated to WRF
#  nomenclature
analysisdat = Dataset(analysispath, 'w')
analysisdat.createDimension('time', size=None)
analysisdat.createDimension('lat', size=ydim)
analysisdat.createDimension('lon', size=xdim)
analysisdat.createVariable(analysisvar, float, ('lat', 'lon'))
analysisdat[analysisvar][:] = analysis
analysisdat.close()
# Store member values to netCDF4 that emulate SENSvals.nc
#  output from sensvector.f
memvalsdat = Dataset(memvalspath, 'w')
memvalsdat.createDimension('xdim', size=xdim)
memvalsdat.createDimension('ydim', size=ydim)
memvalsdat.createDimension('member', size=ensnum)
memvalsdat.createVariable(memvalvar, float, ('member', 'ydim', 'xdim'))
memvalsdat[memvalvar][:] = mem_sens_vals
memvalsdat.close()

def test_point():
    """
    Test to ensure that the point method of subsetting produces
    the location in the sensitivity field with the highest
    sensitivity magnitude and masks all others.
    """
    mask, maxval = esens_subsetting.point(sens_field)
    assert(sens_field[mask == False] == np.max(np.absolute(sens_field)))
    return

def test_ensSubset():
    """
    Tests to ensure sensitivity-weighted and RMS subsetting techniques
    can produce the correct subset given an arbitrary sensitivity field
    and an ensemble where member 5 will always have lowest errors.
    """
    # Attempt with projection method (method=2)
    submems, sensstrings = esens_subsetting.ensSubset(wrfsensfile=senspath,
                                analysis=analysispath, newensnum=3,
                                memvalsfile=memvalspath, fullensnum=ensnum,
                                sensvars=[sensvar], method=2, sens_thresh=0.0)
    print("Subset members in order of increasing error: {} \nSensitivity variables used: {}".format(submems,
            sensstrings))
    # We know that member 5 has the lowest errors over the whole domain,
    #  so it should also have the lowest sensitivity-weighted errors.
    #  Make sure that's the case here.
    assert(submems[0] == 5)

    # Attempt with RMS technique (method=3) and show the same member is
    #  produced even with a nonzero sensitivity threshold.
    submems, sensstrings = esens_subsetting.ensSubset(wrfsensfile=senspath,
                                    analysis=analysispath, newensnum=3,
                                    memvalsfile=memvalspath, fullensnum=ensnum,
                                    sensvars=[sensvar], method=3, sens_thresh=0.8)
    assert(submems[0] == 5)

    return
