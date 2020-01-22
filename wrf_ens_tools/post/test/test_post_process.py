"""
Code for testing routines found in post_process.py
"""
import numpy as np
from numpy.testing import assert_almost_equal

from wrf_ens_tools.post import nearest_lat_lon_index


def test_nearest_lat_lon_index():
    """test the nearest index function for lat/lon coordinates"""
    lat = np.arange(30, 50, 5)
    lon = np.arange(-110, -90, 5)
    ilatlon = [40., -100.]
    clon, clat = np.meshgrid(lon, lat)
    nearest_idx = nearest_lat_lon_index(ilatlon, clat, clon)
    nearest_lat = clat[nearest_idx]
    nearest_lon = clon[nearest_idx]
    assert_almost_equal(nearest_lat, 40.)
    assert_almost_equal(nearest_lon, -100.)
