#!/bin/bash
#SBATCH --job-name=wrf-partmc-les
#SBATCH --nodes=8
#SBATCH -n 384
#SBATCH --partition=sesempi
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=sf20@illinois.edu

# CONFIGURABLE PARAMETERS -------------------------------------------------------------

# Initial condition and emission profile parameters
scenario='uniform-basecase'
#scenario='fx2fy2'
#scenario='fx1fy0'
#scenario='road-10x'
#scenario='point-source-10x10'
#scenario='point-source-1x1'

# Emission rates for aerosols, gases
EMISS_PATH=/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les/emissions
gas_emiss_data_path=$EMISS_PATH/urbanplume_gas_emiss.csv
#gas_emiss_data_path=$EMISS_PATH/urbanplume_gas_emiss_no-ammonia.csv # SF 4/23/24: Run scenario with no ammonia gas emissions 
aero_emiss_data_path=$EMISS_PATH/urbanplume_aero_emiss.csv

# Initial Conditions for aerosols, gases
IC_PATH=/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les/initial-conditions
aero_ic_data_path=$IC_PATH/urbanplume_aero_ics.json
gas_ic_data_path=$IC_PATH/urbanplume_gas_ics.csv
#aero_ic_data_path=$IC_PATH/urbanplume_aero_ics_no-ammonium.json # SF 4/23/24: Run scenario with no initial ammonium in aerosol
#gas_ic_data_path=$IC_PATH/urbanplume_gas_ics_no-ammonia.csv # SF 4/23/24: Run scenario with no initial ammonia gas conc

# Sounding 
SOUNDING_FILE=input_sounding_base
#SOUNDING_FILE=input_sounding_meanwind

# Path to WRF-LES simulation
SIM_PATH=/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les
cd $SIM_PATH

#ARCHIVE_PATH=/data/nriemer/d/sf20/les_output/wrf-partmc
ARCHIVE_PATH=/data/keeling/a/sf20/e/wrf-partmc-gridded-output  # SF 4/30/24: Add path for gridded output directory on E drive
OUTPUT_PATH=$ARCHIVE_PATH/slurm-$SLURM_JOB_ID

# END CONFIGURABLE PARAMETERS-------------------------------------------------------

mkdir $OUTPUT_PATH
mkdir $OUTPUT_PATH/aero_emit_dists
mkdir $OUTPUT_PATH/ics
mkdir $OUTPUT_PATH/out2

cp $SIM_PATH/partmc-files/aero_data.dat $OUTPUT_PATH/aero_data.dat
cp $SIM_PATH/partmc-files/gas_data.dat $OUTPUT_PATH/gas_data.dat
cp $SIM_PATH/partmc-files/gas_params.csv $OUTPUT_PATH/gas_params.csv
cp $SIM_PATH/ideal.exe $OUTPUT_PATH/ideal.exe
cp $SIM_PATH/namelists-and-soundings/$SOUNDING_FILE $OUTPUT_PATH/input_sounding
cp $SIM_PATH/LANDUSE.TBL $OUTPUT_PATH/LANDUSE.TBL
cp $SIM_PATH/namelist.input $OUTPUT_PATH/namelist.input
cp $SIM_PATH/wrf.exe $OUTPUT_PATH/wrf.exe

cd $OUTPUT_PATH

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "Start : $timestamp"
echo ${run}
source ~/.bashrc
#free the system to have unlimited memory requests
ulimit -s -S unlimited
#make permissions group readable on output files
umask 022
#load libraries that you compiled with
module purge
module load gnu/hdf5-1.10.6-gnu-9.3.0
module load gnu/netcdf4-4.7.4-gnu-9.3.0
module load gnu/openmpi-3.1.6-gnu-9.3.0

export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=COMPATIBLE

# Get simulation configuration (chemical mechanism, domain dimensions)
file_path="namelist.input"
# Loop through each line in the file
while IFS= read -r line; do
    #echo $line
  if [[ $line == *"chem_opt "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) # only look at chem_opt for domain 1
    CHEM_OPT=$(echo "$line" | sed 's/[^0-9]*//g') # remove non-integer values
    echo "Using chemical mechanism: $CHEM_OPT"
  fi
  
  if [[ $line == *"s_we "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) # only look at domain 1
    S_WE=$(echo "$line" | sed 's/[^0-9]*//g') # remove non-integer values
  fi
  if [[ $line == *"e_we "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) # only look at domain 1
    E_WE=$(echo "$line" | sed 's/[^0-9]*//g') # remove non-integer values
  fi

  if [[ $line == *"s_sn "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) # only look at domain 1
    S_SN=$(echo "$line" | sed 's/[^0-9]*//g') # remove non-integer values
  fi
  if [[ $line == *"e_sn "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) # only look at domain 1
    E_SN=$(echo "$line" | sed 's/[^0-9]*//g') # remove non-integer values
  fi
  
  if [[ $line == *"s_vert "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) # only look at domain 1
    S_VERT=$(echo "$line" | sed 's/[^0-9]*//g') # remove non-integer values
  fi
  if [[ $line == *"e_vert "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) # only look at domain 1
    E_VERT=$(echo "$line" | sed 's/[^0-9]*//g') # remove non-integer values
  fi
  if [[ $line == *"do_gridded_output "* ]]; then
    line=$(echo "$line" | cut -d ',' -f 1) 
    GRIDDED_OUTPUT=$(echo "$line" | cut -d '=' -f 2 | cut -d ' ' -f 2)
    if [[ $GRIDDED_OUTPUT == ".true." ]]; then
      GRIDDED_OUTPUT=1.0
    else
      GRIDDED_OUTPUT=0.0
    fi
  fi
done < "$file_path"

S_WE=$(expr $S_WE)
E_WE=$(expr $E_WE)
S_SN=$(expr $S_SN)
E_SN=$(expr $E_SN)
S_VERT=$(expr $S_VERT)
E_VERT=$(expr $E_VERT)

extent_we=$(expr $E_WE - $S_WE)
extent_sn=$(expr $E_SN - $S_SN)
extent_vert=$(expr $E_VERT - $S_VERT)

echo "Number of grid cells in west-east: $extent_we"
echo "Number of grid cells in south-north: $extent_sn"
echo "Number of grid cells in vertical: $extent_vert"

# set emissions (gases and aerosols)
python $SIM_PATH/create_aero_emit_dists.py $OUTPUT_PATH $scenario $extent_we $extent_sn $extent_vert $gas_emiss_data_path $aero_emiss_data_path

# set intial conditions (uniform distribution)
python $SIM_PATH/create_aero_ics.py $OUTPUT_PATH 'uniform-basecase' $extent_we $extent_sn $extent_vert $aero_ic_data_path

echo
time mpirun -np 8 ./ideal.exe 

# set gas initial conditions
python $SIM_PATH/write_gas_ic_json.py $OUTPUT_PATH $CHEM_OPT $scenario $gas_ic_data_path
python $SIM_PATH/edit_wrfinput_initcond.py $OUTPUT_PATH $CHEM_OPT $extent_we $extent_sn $extent_vert $scenario

echo 
time mpirun -np $SLURM_NPROCS ./wrf.exe

echo
timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "End : $timestamp"

# move gridded output
if [ $GRIDDED_OUTPUT == 1.0 ]; then
  echo "Moving gridded output files to $OUTPUT_PATH"
  mv $ARCHIVE_PATH/gridded-output*.nc $OUTPUT_PATH
fi

cp $SIM_PATH/slurm-$SLURM_JOB_ID.out $OUTPUT_PATH/slurm-$SLURM_JOB_ID.out
rm $SIM_PATH/slurm-$SLURM_JOB_ID.out

