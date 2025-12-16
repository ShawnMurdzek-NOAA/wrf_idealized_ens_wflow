# Workflow for Idealized WRF Ensembles

Various scripts and utilities for running idealized WRF ensembles.

## Instructions to run

### Setup: Only need to do this once

*In the following sections, {machine} is the machine name. For UNL HCC, {machine} = swan*

1. If it does not already exist for the desired machine, create a WRF environment file in `./env/wrf_{machine}.env`.  
2. Compile WRF using the environment in `./env/wrf_{machine}.env` using the `em_quarter_ss` option (e.g., using `./compile em_quarter_ss`)
3. Create a Python environment using the environment in `./env/py_wrf_ens.yml`. If conda is available, this can be done using `conda env create -f ./env/py_wrf_ens.yml`
4. Create a Python environment file using the Python environment created in (3) in `./env/python_{machine}.env`.

### Running an ensemble

1. If desired, edit the WRF namelist in `./fix/namelist.input.TEMPLATE`. This namelist will be used by all ensemble members. Do not edit anything in `{ }`. Those parameters are placeholders that are changed using `sed` in `run_wrf.submit`.
2. Edit the "User-defined parameters" portion of `run_wrf.submit`. The comments in that script explain the various variables that need to be set.
3. `sbatch run_wrf.submit`

## Brief Description of `run_wrf.submit`

The `run_wrf.submit` script does the following for each ensemble member in serial:

1. Create a run directory in `${run_dir}/${ens_subdir[i]}/run`
2. Modify the WRF namelist using the user-specified options (e.g., `${pbl_opt}`, `${sfclay_opt}`, etc.). This allows the user to run multi-physics and SKEB ensembles.
3. Run `./ideal.exe`
4. If `${mod_wrfinput} -eq 1`, modify `wrfinput_d01` using the user-specified CM1 restart file. This allows the user to change the WRF initial conditions to be output interpolated from CM1.
5. Run `./wrf.exe`

## Resources
- [WRF documentation](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/overview.html)
- [Description of WRF namelist options](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/namelist_variables.html)
