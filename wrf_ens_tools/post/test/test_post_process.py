"""
Code for testing routines found in post_process.py
"""
import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from wrf_ens_tools.post import destagger, earth_relative_winds, nearest_lat_lon_index, simulated_reflectivity


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


def test_destagger():
    """test the destagger function"""
    flat = np.linspace(0, 20, 12)
    u_stag = flat.reshape((-1, 3, 2))
    u = destagger(u_stag, 1)
    u_ref = np.array([[[1.81818182, 3.63636364],
                       [5.45454545, 7.27272727]],
                      [[12.72727273, 14.54545455],
                       [16.36363636, 18.18181818]]])
    assert_array_almost_equal(u_ref, u)


def test_earth_relative_winds():
    """test for earth relative wind computation with no rotation"""
    u = 10.
    v = 10.
    sinalpha = 0.
    cosalpha = 1.
    u_rot, v_rot = earth_relative_winds(u, v, sinalpha, cosalpha)
    assert_almost_equal(u_rot, u)
    assert_almost_equal(v_rot, v)


def test_earth_relative_winds_rotation():
    """test for earth relative wind computation with no rotation"""
    u = 10.
    v = 10.
    sinalpha = 0.015
    cosalpha = 0.995
    u_rot, v_rot = earth_relative_winds(u, v, sinalpha, cosalpha)
    assert_almost_equal(u_rot, 9.799999999999999)
    assert_almost_equal(v_rot, 10.1)


def test_dbz():
    """test the simulated reflectivity function"""
    t = np.linspace(250, 300, 4).reshape((2, 2))
    p = np.linspace(99000, 101000, 4).reshape((2, 2))
    q = np.linspace(0.003, 0.005, 4).reshape((2, 2))
    dbz = simulated_reflectivity(p, t, q, q)
    truth = np.array([[53.88258, 54.96516],
                     [55.821661, 56.522127]])
    assert_array_almost_equal(dbz, truth)