"""
Convert UAS CSV Files to DART obs_seq.out File Format

Assume the following units for the inputs CSV file:
    t: seconds since CM1 NR started (1200 UTC 15 April 2009)
    x, y, z: km, with the origin being the center of the CM1 domain
    P: Unused, so units are not important
    T: K
    U, V: m/s
    W: Unused, so units are not important

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import datetime as dt
import sys
import argparse
import pyproj


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

    parser = argparse.ArgumentParser(description='Convert a UAS CSV observation file to an \
                                                  obs_seq.out file that can be used by DART.')
    
    # Positional arguments
    parser.add_argument('in_file', 
                        help='Input UAS CSV file',
                        type=str)
    
    # Optional arguments
    parser.add_argument('-p',
                        dest='proj_str',
                        default='+proj=lcc +lat_0=40 +lon_0=-100 +lat_1=40 +lat_2=40',
                        help='Map projection string. See \
                              https://proj.org/en/stable/operations/projections/lcc.html for details',
                        type=str)
    
    parser.add_argument('-t',
                        dest='ref_time',
                        default='20090415120000',
                        help='Reference time. All times in the CSV file (column "t") are seconds \
                              since the reference time. Format: YYYYMMDDHHMMSS',
                        type=str)
    
    return parser.parse_args(argv)


def read_uas_csv(fname, 
                 ref_time=dt.datetime(2009, 4, 15, 12),
                 proj_str='+proj=lcc +lat_0=40 +lon_0=-100 +lat_1=40 +lat_2=40'):
    """
    Read UAS CSV file and perform necessary conversions

    Parameters
    ----------
    fname : string
        UAS CSV file to read
    ref_time : dt.datetime, optional
        Reference time used by the CSV file. All times in the CSV are seconds since the ref_time.
    proj_str : string
        Map projection used to convert UAS (x, y) coordinates to (lat, lon) coordinates
        See https://proj.org/en/stable/operations/projections/lcc.html for details.

    Returns
    -------
    out_df : pd.DataFrame
        CSV file output

    """
    
    # Read CSV
    out_df = pd.read_csv(fname, skiprows=1)
    
    # Convert time to days and seconds since 0000 UTC 1 Jan 1601
    start_t = (ref_time - dt.datetime(1601, 1, 1, 0)).total_seconds()
    total_t = start_t + out_df['t'].values
    s_in_day = 24*60*60
    out_df['dart_days'] = np.int64(np.floor(total_t / s_in_day))
    out_df['dart_sec'] = np.int64(total_t % s_in_day)
    
    # Convert (x, y) to (lat, lon) coordinates using map projection
    proj = pyproj.Proj(param.proj_str)
    lon, lat = proj(out_df['x'], out_df['y'], inverse=True)
    lon[lon < 0] = lon[lon < 0] + 360.
    out_df['dart_lon'] = np.deg2rad(lon)
    out_df['dart_lat'] = np.deg2rad(lat)
    
    # Convert height from km to m
    out_df['dart_z'] = out_df['z'] * 1000.
    
    return out_df


def write_uas_df_to_obs_seq(df, obserr={'T':1, 'U':1, 'V':1}):
    """
    Write a UAS obs DataFrame to a DART obs_seq.out file

    Parameters
    ----------
    df : pd.DataFrame
        UAS observations
    obserr : dictionary
        Observation error variance for each observation column in the DataFrame

    Returns
    -------
    None.

    """
    
    # Observation types
    # Key: Column name in DataFrame
    # name: Variable name in DART
    # type: Variable type in DART
    ob_types = {'T': {'name': 'UAS_TEMPERATURE', 'type': 118},
                'U': {'name': 'UAS_U_WIND_COMPONENT', 'type': 116},
                'V': {'name': 'UAS_V_WIND_COMPONENT', 'type': 117}}
    
    # Determine the number of observations
    # NaN can be used for missing observations
    nobs = 3*len(df)
    for v in ob_types:
        nobs = nobs - np.sum(np.isnan(df[v]))
    
    # Write obs_seq.out file
    fname = 'obs_seq.out'
    with open(fname, 'w') as fptr:
        fptr.write(" obs_sequence\n")
        fptr.write("obs_kind_definitions\n")

        fptr.write(f"    {len(ob_types)} \n")
        for v in ob_types:
            fptr.write(f"    {ob_types[v]['type']}          {ob_types[v]['name']}   \n")
    
        fptr.write("  num_copies:            1  num_qc:            1\n")
        fptr.write(f" num_obs:       {nobs}  max_num_obs:       {nobs}\n")
        fptr.write("MADIS observation\n")
        fptr.write("Data QC\n")
        fptr.write(f"  first:            1  last:       {nobs}\n")
        
        t_sec = df['dart_sec'].values
        t_day = df['dart_days'].values
        lon = df['dart_lon'].values
        lat = df['dart_lat'].values
        z = df['dart_z'].values
        n = 0
        truth = 1
        for v in ob_types:
            kind = ob_types[v]['type']
            oe = obserr[v]
            obs = df[v].values
            for i, val in enumerate(obs):
                fptr.write(f" OBS            {n+1}\n")
                fptr.write(f"   {val:20.14f}\n")
                fptr.write(f"   {truth:20.14f}\n")

                if nobs == 1:
                    fptr.write(" -1 -1 -1\n") # Only 1 ob
                elif n+1 == 1: 
                    fptr.write(f" -1 {n+2} -1\n") # First ob
                elif n+1 == nobs:
                    fptr.write(f" {n} -1 -1\n") # Last ob
                else:
                    fptr.write(f" {n} {n+2} -1\n") 
            
                fptr.write("obdef\n")
                fptr.write("loc3d\n")
            
                fptr.write(f"    {lon[i]:20.14f}          {lat[i]:20.14f}          {z[i]:20.14f}     3\n")
                fptr.write("kind\n")       
                fptr.write(f"     {kind}     \n")       
                fptr.write(f"    {t_sec[i]}          {t_day[i]}     \n")               
                fptr.write(f"    {oe:20.14f}  \n")
                
                n = n + 1
                
    return None  


if __name__ == '__main__':
    
    start = dt.datetime.now()
    print('Starting convert_csv_to_obs_seq.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Read in parameters and UAS CSV file
    param = parse_in_args(sys.argv[1:])
    ref_time = dt.datetime.strptime(param.ref_time, '%Y%m%d%H%M%S')
    uas_df = read_uas_csv(param.in_file, ref_time=ref_time, proj_str=param.proj_str)
    
    # Write DataFrame to DART obs_seq.out file
    write_uas_df_to_obs_seq(uas_df, obserr={'T':1, 'U':1, 'V':1})

    print('Program finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End convert_csv_to_obs_seq.py
"""