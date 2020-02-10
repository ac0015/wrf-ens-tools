import numpy as np
import xarray as xr
from wrf_ens_tools.calc import *
from wrf_ens_tools.post import storeIdealizedPracPef, storeIdealizedNearestNeighborFortran

# Create dummy data
sixhr = True; timedim = 6
nbrhd = 12.; dx = 4.; sigma = 1
xdim, ydim, zdim = 588, 540, 40
ypts, xpts = np.meshgrid(np.arange(ydim), np.arange(xdim))
fake_uh_fcst = np.random.uniform(low=0., high=100., size=(timedim, ydim, xdim))
fake_refl_fcst = np.zeros((timedim, zdim, ydim, xdim))

# Fill one vertical column with 45 dBZ reflectivity values at arbitrary time slice
fake_refl_fcst[int(timedim/2), :, int(ydim/2), int(xdim/2)] = 45.
# Number of grid pts around uh pts where refl exceedance is checked
#  (see Sobash et al 2011)
r = 25. / dx

# Mask used to determine number of total pts within refl neighborhood radius
mask = dist_mask(xind=int(xdim/2), yind=int(ydim/2),
                xpts=xpts, ypts=ypts, r=25.)
max_num_ssr = len(mask == True)
fhr = 28 # random forecast hour to attribute to meta data

# Get lat/lons for interpolation from SPC 211 grid
ds = xr.open_dataset("wrfoutREFd2")
lats = ds["XLAT"][0]
lons = ds["XLONG"][0]
ds.close()

def test_ssr_sspf_spc_grid():
    # Create SSRs and SSPFs from dummy data
    # SSR/SSPF arrays on SPC 211 grid
    ssr211 = gen_surrogate_severe_reports(uh_arr=fake_uh_fcst,
                                    sim_refl_arr=fake_refl_fcst,
                                    uh_thresh=40., lats=lats,
                                    lons=lons, dx=dx, spc_grid=True)
    sspf211 = gen_SSPFs_from_SSRs(ssr_arr=ssr211, sigma=sigma)
    # Store SSPF and SSR data to netCDF files
    storeIdealizedPracPef(sspf_arr=sspf211, outlats=lats, outlons=lons,
                            outpath="idealized_pperf_from_spc211grid.nc",
                            sigma=sigma, fhrs=fhr, spc_grid=True)
    pp_ds = xr.open_dataset("idealized_pperf_from_spc211grid.nc")
    pperf211 = pp_ds["practically_perfect"][:]
    num_ssrs = len(np.isclose(ssr211, 1) == True)
    # Ensure SSRs are generating correctly
    assert(num_ssrs <= max_num_ssr)
    # # Ensure SSPFs are valid
    # assert(pperf211.all() < 0.99)

def test_ssr_sspf_native_grid():
    # Create SSRs and SSPFs from dummy data
    # SSR/SSPF arrays on native WRF grid
    ssr = gen_surrogate_severe_reports(uh_arr=fake_uh_fcst,
                                sim_refl_arr=fake_refl_fcst,
                                uh_thresh=40., lats=lats, lons=lons,
                                dx=dx, spc_grid=False)
    sspf = gen_SSPFs_from_SSRs(ssr_arr=ssr, sigma=sigma)
    # Store SSPF and SSR data to netCDF files
    storeIdealizedPracPef(sspf_arr=sspf, outlats=lats, outlons=lons,
                            outpath="idealized_pperf_from_nativegrid.nc",
                            sigma=sigma, fhrs=fhr, spc_grid=False)
    storeIdealizedNearestNeighborFortran(ssr_arr=ssr,
                            outpath="idealized_rel_ob_grid.nc",
                            wrfrefpath="wrfoutREFd2",
                            obvar="P_HYD")
    num_ssrs = len(np.isclose(ssr, 1) == True)
    stored_ssrs = xr.open_dataset("idealized_rel_ob_grid.nc")
    num_stored_ssrs = len(np.isclose(stored_ssrs["P_HYD"][0,0]) == True)
    # Ensure reliability ob storage is working correctly
    assert(num_ssrs <= max_num_ssr)
    # assert(ssr == stored_ssrs["P_HYD"][0,0])
    # stored_ssrs.close()
    # # Ensure SSPF is working correctly
    # pp_ds = xr.open_dataset("idealized_pperf_from_nativegrid.nc")
    # pperf = pp_ds["practically_perfect"][:]
    # pp_ds.close()
    # assert(sspf.all() == pperf.all())
    # assert(sspf < 0.99)
