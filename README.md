# Workflow for Idealized WRF Ensembles

Various scripts and utilities for running idealized WRF ensembles.

## Instructions to run

### Setup: Only need to do this once

*In the following sections, {machine} is the machine name. For UNL HCC, {machine} = swan*

1. If it does not already exist for the desired machine, create a WRF environment file in `./env/wrf_{machine}.env`.  
2. Compile WRF using the environment in `./env/wrf_{machine}.env`.
3. Create a Python environment using the environment in `./env/py_wrf_ens.yml`. If conda is available, this can be done using `conda env create -f ./env/py_wrf_ens.yml`
4. Create a Python environment file using the Python environment created in (3) in `./env/python_{machine}.env`.

### Running an ensemble

1. Edit the "User-defined parameters" portion of `run_wrf.submit`
2. `sbatch run_wrf.submit`

## Description of User-Defined Parameters

### Paths
- `wrf_env`: File containing environment used to run WRF.
- `wrf_src`: Directory containing compiled idealized (em_quarter_ss) WRF code.
- `run_dir`: Parent directory where are ensemble runs will be performed.
- `cm1_rst`: CM1 restart file used to initialize the idealized WRF runs.

### Ensemble member configurations
*Each variable is a list with each entry corresponding to a different ensemble member*  
- `ens_subdir`: Name of the subdirectory for each ensemble member.  
- `pbl_opt`: PBL scheme option.  
- `sfclay_opt`: Surface layer scheme option.  
- `ra_sw_opt`: Shortwave radiation scheme option.  
- `ra_lw_opt`: Longwave radiation scheme option.  
- `skebs_opt`: Option to use SKEBs.  

### Other
- `mod_wrfinput`: Option to modify `wrfinput_d01` using the `cm1_rst` restart file  
- `ndcnst`: Constant droplet number concentration (cm^-3) used by py_scripts/wrf/modify_wrf_with_cm1.py  

## Resources
- [WRF documentation](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/overview.html)
