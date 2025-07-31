# Workflow for Idealized WRF Ensembles on UNL HCC

## Instructions to Run

1. Compile WRF
2. Edit the "User-defined parameters" portion of `run_wrf.submit`
3. `sbatch run_wrf.submit`

## Description of User-Defined Parameters

### Paths
- `wrf_env`: File containing environment used to run WRF.
- `wrf_src`: Directory containing compiled idealized (em_quarter_ss) WRF code.
- `run_dir`: Parent directory where are ensemble runs will be performed.
- `cm1_rst`: CM1 restart file used to initialize the idealized WRF runs.
- `py_script_dir`: Path to [py_scripts](https://github.com/ShawnMurdzek-NOAA/py_scripts/tree/main). Needed to modify wrfinput file with CM1 restart fields.

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
