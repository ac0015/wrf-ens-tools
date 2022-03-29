from setuptools import setup

setup(
    name='wrf-ens-tools',
    version='0.1',
    install_requires=['netcdf4', 'numpy', 'scipy', 'xarray', 'pandas',
                      'wrf-python', 'matplotlib', 'siphon', 'cartopy', 'sklearn',
                      'dask', 'metpy>=1.0', 'cmocean'],
    license='BSD-3',
    author='Austin Coleman, Russel Manser, Tyler Wixtrom',
)
