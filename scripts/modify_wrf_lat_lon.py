"""
Modify WRF Latitude and Longitude Coordinates Using a Map Projection

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import sys
import argparse
import xarray as xr
import pyproj
import numpy as np


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def parse_in_args(argv):
    """
    Parse input arguments

    Parameters
    ----------
    argv : list
        Command-line arguments from sys.argv[1:]
    
    Returns
    -------
    Parsed input arguments

    """

    parser = argparse.ArgumentParser(description='Modify WRF (lat, lon) coordinates. Input file \
                                                  will be modified in place (i.e., the XLAT and \
                                                  XLONG fields will be overwritten)')
    
    # Positional arguments
    parser.add_argument('in_file', 
                        help='Input WRF netCDF file (e.g., a wrfout file)',
                        type=str)

    # Optional arguments
    parser.add_argument('-p',
                        dest='proj_str',
                        default='+proj=lcc +lat_0=40 +lon_0=-100 +lat_1=40 +lat_2=40',
                        help='Map projection string. See \
                              https://proj.org/en/stable/operations/projections/lcc.html for details',
                        type=str)
    
    return parser.parse_args(argv)



def compute_wrf_grid(wrf_ds):
    """
    Create an (x, y, z) grid for WRF output in m
    """

    wrf_grid = {}

    # X direction
    wrf_grid['west_east'] = wrf_ds.attrs['DX'] * wrf_ds['west_east'].values
    wrf_grid['west_east_stag'] = wrf_ds.attrs['DX'] * wrf_ds['west_east_stag'].values
    wrf_grid['west_east_stag'] = wrf_grid['west_east_stag'] - 0.5*(wrf_grid['west_east'][-1] + wrf_ds.attrs['DX'])
    wrf_grid['west_east'] = wrf_grid['west_east'] - 0.5*wrf_grid['west_east'][-1]

    # Y direction
    wrf_grid['south_north'] = wrf_ds.attrs['DY'] * wrf_ds['south_north'].values
    wrf_grid['south_north_stag'] = wrf_ds.attrs['DY'] * wrf_ds['south_north_stag'].values
    wrf_grid['south_north_stag'] = wrf_grid['south_north_stag'] - 0.5*(wrf_grid['south_north'][-1] + wrf_ds.attrs['DY'])
    wrf_grid['south_north'] = wrf_grid['south_north'] - 0.5*wrf_grid['south_north'][-1]

    return wrf_grid


def compute_lat_lon_coords(wrf_grid, proj):
    """
    Compute (lat, lon) coordinates for each Cartesian coordinate using the provided map projection
    """
    
    wrf_lat_lon = {}
    
    lat_lon_names = [['XLONG', 'XLAT'],
                     ['XLONG_U', 'XLAT_U'],
                     ['XLONG_V', 'XLAT_V']]
    cart_names = [['west_east', 'south_north'],
                  ['west_east_stag', 'south_north'],
                  ['west_east', 'south_north_stag']]
    
    for n1, n2 in zip(lat_lon_names, cart_names):
        x, y = np.meshgrid(wrf_grid[n2[0]], wrf_grid[n2[1]])
        wrf_lat_lon[n1[0]], wrf_lat_lon[n1[1]] = proj(x, y, inverse=True)

    return wrf_lat_lon


def overwrite_wrf_lat_lon(wrf_ds, wrf_lat_lon):
    """
    Overwrite the latitude and longitude fields in the input DataSet
    """
    
    fields = ['XLAT', 'XLONG', 
              'XLAT_U', 'XLONG_U',
              'XLAT_V', 'XLONG_V']
    for f in fields:
        wrf_ds[f].values[0, :, :] = wrf_lat_lon[f]
    
    return wrf_ds


if __name__ == '__main__':
    
    start = dt.datetime.now()
    print('Starting modify_wrf_lat_lon.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Read in parameters, netCDF file, and create projection
    param = parse_in_args(sys.argv[1:])
    ds = xr.open_dataset(param.in_file, mode='a')
    proj = pyproj.Proj(param.proj_str)
    
    # Create WRF grid in meters
    cart_grid = compute_wrf_grid(ds)
    
    # Compute updated XLAT and XLONG fields
    sphere_grid = compute_lat_lon_coords(cart_grid, proj)
    ds = overwrite_wrf_lat_lon(ds, sphere_grid)
    
    # Save output
    ds.to_netcdf(param.in_file, mode='a')

    print('Program finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End modify_wrf_lat_lon.py
"""