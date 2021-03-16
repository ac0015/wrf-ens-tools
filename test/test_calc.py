# test_calc.py

"""Test the calc module"""

import xarray as xr
from numpy.testing import assert_almost_equal

from wrf_ens_tools.calc import FSS

def test_FSS():
    """Test the FSS function."""
    fcst = xr.open_dataset('test/forecast_probs.nc')
    obs = xr.open_dataset('test/observation_probs.nc')

    actual_fss, actual_fbs, actual_fbs_ref = FSS(
        fcst.REFL_10CM_25.values,
        obs.col_max_refl_25.values,
        return_fbs=True
    )

    true_fss = 0.59736522
    true_fbs = 608.62545956
    true_fbs_ref = 1511.60675737

    assert_almost_equal(actual_fss, true_fss)
    assert_almost_equal(actual_fbs, true_fbs)
    assert_almost_equal(actual_fbs_ref, true_fbs_ref)
