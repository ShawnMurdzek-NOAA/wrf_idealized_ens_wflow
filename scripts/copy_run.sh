
# Script to create a clean run directory

# Passed arguments:
#     $1 : WRF source code directory
#     $2 : Destination for run directory

wrf_src=$1
run_dir=$2

echo
echo "Inside copy_run.sh"
echo "=================="
echo "WRF directory = ${wrf_src}"
echo "Run directory = ${run_dir}"

if [[ -d ${run_dir}/run ]]; then
  echo
  echo 'Deleting old run directory'
  echo
  rm -r ${run_dir}/run
fi

cp -r ${wrf_src}/run ${run_dir}/
cd ${run_dir}/run

# Remove dangling links
rm ideal.exe
rm input_sounding
rm MPTABLE.TBL
rm wrf.exe

# Copy over executables
cp ${wrf_src}/main/ideal.exe .
cp ${wrf_src}/main/wrf.exe .
cp ${wrf_src}/phys/noahmp/parameters/MPTABLE.TBL .

echo "copy_run.sh done"
echo
