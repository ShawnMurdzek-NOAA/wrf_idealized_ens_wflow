"""
Quick script to add global attributes to a netCDF file

Passed Arguments
----------------
sys.argv[1] : string
    netCDF file name

shawn.s.murdzek@noaa.gov
"""

import xarray as xr
import sys

attrs = {'CEN_LAT': 40.,
         'CEN_LON': -100.,
         'TRUELAT1': 40.,
         'TRUELAT2': 40.,
         'MOAD_CEN_LAT': 40.,
         'STAND_LON': -100.,
         'POLE_LAT': 90.,
         'POLE_LON': 0.,
         'MAP_PROJ': 1,
         'MAP_PROJ_CHAR': "Lambert Conformal"}

nc_fname = sys.argv[1]

ds = xr.open_dataset(nc_fname, mode='a')

for key in attrs:
    ds.attrs[key] = attrs[key]

ds.to_netcdf(nc_fname, mode='a')
