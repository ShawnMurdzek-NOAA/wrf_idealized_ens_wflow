
#---------------------------------------------------------------------------------------------------
# Description

# This initializes various environment variables used by the slurm submission scripts (*.submit)

#---------------------------------------------------------------------------------------------------

# Paths
# -----
# run_dir: Where the WRF ensemble will be run
# wrf_src: Location of the compiled WRF code
# cm1_rst: CM1 restart file used to modify wrfinput_d01. Only used if mod_wrfinput=1
# home: Should not need to change this
export RUN_DIR=/work/ahouston/smurdzek/wrf_ens_github_test
export WRF_SRC=/home/ahouston/smurdzek/WRF_code/idealized/em_quarter_ss_modified/WRF
export CM1_RST=/work/ahouston/smurdzek/CM1_output/cm1out_rst_000004.nc
export WFLOW_HOME=`pwd`

# Machine (needed for loading the proper environment)
# ---------------------------------------------------
# Note that the top portion of the slurm submission scripts (i.e., the SBATCH parameters) will 
# likely need to be edited for any machine that is not "swan"
export MACHINE='swan'
