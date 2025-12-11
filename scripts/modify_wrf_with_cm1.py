"""
Modify WRF Input File Using CM1 Output

shawn.murdzek@colorado.edu
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import sys
import argparse
import copy
import numpy as np
import xarray as xr
import pandas as pd
import metpy.calc as mc
from metpy.units import units
import scipy.interpolate as si


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

    parser = argparse.ArgumentParser(description='Script that interpolates CM1 fields to a \
                                                  wrfinput file')
    
    # Positional arguments
    parser.add_argument('cm1_file', 
                        help='CM1 netCDF restart file (e.g., cm1out_rst_XXXXXX.nc)',
                        type=str)
    
    parser.add_argument('cm1_base', 
                        help='CM1 input_sounding text file',
                        type=str)

    parser.add_argument('wrf_infile', 
                        help='Input wrfinput netCDF file',
                        type=str)

    parser.add_argument('wrf_outfile', 
                        help='Output wrfinput netCDF file',
                        type=str)

    # Optional arguments
    parser.add_argument('-i',
                        dest='iopt',
                        default='nearest',
                        help='Interpolation option. Passed to scipy.RegularGridInterpolator',
                        type=str)

    parser.add_argument('--ndcnst',
                        dest='ndcnst',
                        default=300.,
                        help='Constant droplet number concentration (cm^-3) used in CM1 \
                              (e.g., for Morrison 2-moment scheme). Set to 0 if droplet number \
                              concentration is predicted in CM1',
                        type=float)

    parser.add_argument('--use_thm',
                        dest='use_thm',
                        default=1,
                        help='Option to use moist theta (THM) as the prognostic temperature \
                              variable in WRF. Value should match use_theta_m in WRF namelist, \
                              which is 1 (True) by default in WRFv4. To not use THM, set to 0.',
                        type=int)

    return parser.parse_args(argv)


def constants():
    """
    Define constants as global variables
    """

    # Taken from WRF/share/module_model_constants.F
    global r_d, cp, r_v, cv, cvpm, p1000mb, t0, rvovrd
    r_d = 287.
    cp = 7. * r_d / 2
    r_v = 461.6
    cv = cp - r_d
    cvpm = -cv / cp
    p1000mb = 100000.
    t0 = 300.
    rvovrd = r_v / r_d

    return None


def read_cm1_restart(fname):
    """
    Read in CM1 restart file
    """

    cm1_ds = xr.open_dataset(fname, decode_timedelta=False)

    # Remove extra times
    cm1_ds = cm1_ds.sel(time=slice(cm1_ds['time'].values[0]))

    return cm1_ds


def read_cm1_input_sounding(fname):
    """
    Read in CM1 base state from input_sounding
    """

    return pd.read_csv(fname, 
                       sep='\s+', 
                       names=['hgt', 'theta', 'qv', 'u', 'v'],
                       skiprows=1)


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
    
    # Z direction
    # Assume geopotential is constant at each vertical level (i.e., there is no terrain)
    geopotential = wrf_ds['PH'].values[0, :, 0, 0] + wrf_ds['PHB'].values[0, :, 0, 0]
    hgt = mc.geopotential_to_height(geopotential * units.m**2 / units.s**2).to(units.m).magnitude
    wrf_grid['bottom_top_stag'] = hgt
    wrf_grid['bottom_top'] = 0.5 * (wrf_grid['bottom_top_stag'][:-1] + wrf_grid['bottom_top_stag'][1:])

    return wrf_grid


def interp_1_field(wrf_ds, cm1_ds, wrf_field, cm1_field, wrf_grid, method='nearest'):
    """
    Interpolate a single field from CM1 to WRF
    """

    # Create RegularGridInterpolator using Scipy
    # Need 0 index to eliminate time dimension (which is always first)
    cm1_coords = []
    noextrapolate = True
    for d in cm1_ds[cm1_field].dims:
        if d == 'time': 
            continue
        elif d == 'ni':
            key = 'xh'
        elif d == 'nip1':
            key = 'xf'
        elif d == 'nj':
            key = 'yh'
        elif d == 'njp1':
            key = 'yf'
        elif d == 'nk':
            key = 'zh'
        elif d == 'nkp1':
            key = 'zf'
            noextrapolate = False  # Vertical extrapolation is often necessary
        cm1_coords.append(cm1_ds[key].values)
    interp = si.RegularGridInterpolator(cm1_coords, cm1_ds[cm1_field][0].values, method=method,
                                        fill_value=None, bounds_error=noextrapolate)

    # Interpolate to WRF grid
    wrf_dims = wrf_ds[wrf_field].dims
    if len(wrf_dims) == 3:
        wrf_coords = np.meshgrid(wrf_grid[wrf_dims[1]], wrf_grid[wrf_dims[2]], indexing='ij')
    elif len(wrf_dims) == 4:
        wrf_coords = np.meshgrid(wrf_grid[wrf_dims[1]], wrf_grid[wrf_dims[2]], wrf_grid[wrf_dims[3]], indexing='ij')
    else:
        raise ValueError(f"WRF field with {len(wrf_dims) - 1} dimensions is not supported")
    out = interp(wrf_coords)
    wrf_ds[wrf_field].values = out[np.newaxis, :]

    return wrf_ds


def interp_cm1_to_wrf(wrf_ds, cm1_ds, cm1_base, wrf_grid, fields, method='nearest', verbose=1):
    """
    Interpolate several fields from CM1 to WRF
    """

    for f in fields.keys():

        if verbose > 0: print(f"Interpolating {f}")

        # Preprocess
        if fields[f] == 'tha':
            # Convert perturbation potential temperature to regular potential temperature
            cm1_ds['theta'] = copy.deepcopy(cm1_ds['tha'])
            th_base = np.interp(cm1_ds['zh'].values, cm1_base['hgt'].values, cm1_base['theta'].values)
            cm1_ds['theta'].values = cm1_ds['tha'].values + th_base[np.newaxis, :, np.newaxis, np.newaxis]
            cm1_ds['theta'].attrs['long_name'] = 'potential temperature'
            fields[f] = 'theta'
        elif fields[f] == 'thm':
            # Compute perturbation moist potential temperature
            cm1_ds['thm'] = copy.deepcopy(cm1_ds['tha'])
            th_base = np.interp(cm1_ds['zh'].values, cm1_base['hgt'].values, cm1_base['theta'].values)
            theta = cm1_ds['tha'].values + th_base[np.newaxis, :, np.newaxis, np.newaxis]
            cm1_ds['thm'].values = (1. + (r_v / r_d) * cm1_ds['qv'].values) * theta - t0
            cm1_ds['thm'].attrs['long_name'] = 'perturbation moist potential temperature'

        # Interpolate
        wrf_ds = interp_1_field(wrf_ds, cm1_ds, f, fields[f], wrf_grid, method=method)

        # Postprocess (e.g., convert from "full" value back to perturbation)
        if f == 'T':
            # Convert regular potential temperature to perturbation potential temperature
            # Value of 300K is hard coded in share/module_model_constants.F in WRF
            wrf_ds['T'].values = wrf_ds['T'].values - t0

    return wrf_ds


def interp_const_droplet_number(wrf_ds, ndcnst):
    """
    Interpolate constant droplet number concentration from CM1 to WRF

    Must interpolate cloud water mixing ratio to WRF first
    """

    # Initialize droplet number mixing ratio to 0
    wrf_ds['QNDROP'].values = np.zeros(wrf_ds['QNDROP'].shape)
    
    # Determine cloudy gridpoints
    icloud = wrf_ds['QCLOUD'].values > 0

    # Compute QNDROP from ndcnst
    wrf_ds['QNDROP'].values[icloud] = ndcnst * 1e6 * (wrf_ds['AL'].values[icloud] + wrf_ds['ALB'].values[icloud])

    return wrf_ds


def recompute_wrf_density(wrf_ds):
    """
    Recompute WRF inverse density perturbations

    Uses the same code as WRF/dyn_em/module_initialize_ideal.F (lines 1172-1217)
    """

    # Compute qv correction for virtual temperature
    qvf = 1. + (rvovrd * wrf_ds['QVAPOR'].values)

    # Compute total inverse density
    alt = ((r_d/p1000mb) * (wrf_ds['T'].values + t0) * qvf * 
           ((wrf_ds['P'].values + wrf_ds['PB'].values) / p1000mb)**cvpm)

    # Compute perturbation inverse density
    wrf_ds['AL'].values = alt - wrf_ds['ALB'].values

    return wrf_ds


def rebalance_PH_hydrostatic(wrf_ds):
    """
    Rebalance perturbation geopotential. Must recompute inverse density perturbations first. 

    Uses the same code as WRF/dyn_em/module_initialize_ideal.F (lines 1172-1217)
    """

    for k in range(1, wrf_ds.attrs['BOTTOM-TOP_GRID_DIMENSION']):
        wrf_ds['PH'].values[:, k, :, :] = wrf_ds['PH'].values[:, k-1, :, :] - (wrf_ds['DNW'].values[:, k-1] *
                                           (((wrf_ds['C1H'].values[:, k-1] * wrf_ds['MUB'].values + wrf_ds['C2H'].values[:, k-1]) +
                                             wrf_ds['C1H'].values[:, k-1] * wrf_ds['MU'].values) *
                                            wrf_ds['AL'].values[:, k-1, :, :] + 
                                            (wrf_ds['C1H'].values[:, k-1] * wrf_ds['MU'].values) * wrf_ds['ALB'].values[:, k-1, :, :]))

    return wrf_ds


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting modify_wrfinput_with_cm1.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Read in data and compute WRF grid
    constants()
    param = parse_in_args(sys.argv[1:])
    cm1_ds = read_cm1_restart(param.cm1_file)
    wrf_ds = xr.open_dataset(param.wrf_infile)
    cm1_base = read_cm1_input_sounding(param.cm1_base)
    wrf_grid = compute_wrf_grid(wrf_ds)

    # Fields to interpolate
    # Note that we need to interpolate T for rebalancing
    fields = {'U':'ua',
              'V':'va',
              'W':'wa',
              'T':'tha',
              'QVAPOR':'qv',
              'QCLOUD':'qc',
              'QRAIN':'qr',
              'QICE':'qi',
              'QSNOW':'qs',
              'QGRAUP':'qg',
              'QNRAIN':'ncr',
              'QNICE':'nci',
              'QNSNOW':'ncs',
              'QNGRAUPEL':'ncg'}

    # Add THM if use_thm = True
    if (param.use_thm == 1):
        fields['THM'] = 'thm'

    # Add cloud droplet number mixing ratio
    if np.isclose(param.ndcnst, 0):
        fields['QNDROP'] = 'nc'

    # Interpolate
    wrf_ds = interp_cm1_to_wrf(wrf_ds, cm1_ds, cm1_base, wrf_grid, fields, method=param.iopt)
    if param.ndcnst > 0:
        print('Interpolating QNDROP')
        wrf_ds = interp_const_droplet_number(wrf_ds, param.ndcnst)

    # Set T = THM if use_thm = False
    if (param.use_thm == 0):
        wrf_ds['THM'].values = wrf_ds['T'].values

    # Recompute inverse perturbation density + rebalance hydrostatically
    print('Recomputing inverse perturbation density and perturbation geopotential')
    wrf_ds = recompute_wrf_density(wrf_ds)
    wrf_ds = rebalance_PH_hydrostatic(wrf_ds)

    # Write out wrf_ds
    wrf_ds.to_netcdf(param.wrf_outfile)

    print('Program finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End modify_wrfinput_with_cm1.py
"""
