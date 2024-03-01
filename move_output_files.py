import os, sys
from datetime import datetime
import shutil
import fnmatch

wrk_dir = os.getcwd()
output_dir = r"/data/nriemer/d/sf20/les_output/wrf-partmc"
now = datetime.now()
datestamp = now.strftime('%Y-%m-%d')
timestamp = datetime.now().strftime('%H%M%S')

slurm_files = fnmatch.filter(os.listdir('.'), "slurm-*.out*")
if len(slurm_files) == 1:
    output_subdir = slurm_files[0].replace('.out', '')
else:
    datestamp = now.strftime('%Y-%m-%d')
    timestamp = datetime.now().strftime('%H%M%S')
    output_path = os.path.join(output_dir, datestamp, timestamp)

output_path = os.path.join(output_dir, output_subdir)
if not os.path.isdir(output_path):
    os.makedirs(output_path)
    print(f'{output_path}')

error_files = fnmatch.filter(os.listdir('.'), "rsl.error.*")
out_files = fnmatch.filter(os.listdir('.'), "rsl.out.*")
wrfinput_files = fnmatch.filter(os.listdir('.'), "wrfinput_*")
wrfoutput_files = fnmatch.filter(os.listdir('.'), "wrfout_*")
wrfrst_files = fnmatch.filter(os.listdir('.'), "wrfrst_*")
aerosol_dist_files = fnmatch.filter(os.listdir('.'), "aerosol_dist_d*")
aerosol_files = fnmatch.filter(os.listdir('.'), "aerosols_d*")

namelist_output = fnmatch.filter(os.listdir('.'), "namelist.output")
species_config_files = fnmatch.filter(os.listdir('.'), "species_*data.json")

move_files = (error_files + out_files + slurm_files + wrfinput_files 
              + wrfoutput_files + wrfrst_files + aerosol_dist_files + 
              aerosol_files + namelist_output + species_config_files)

if len(move_files) == 0:
    print('No files to move, exiting')
    sys.exit()

print('Moving simulation files')
for file in move_files:
    shutil.move(f"{wrk_dir}/{file}", f"{output_path}/{file}")

wrfpartmc_output_path = '/data/nriemer/d/sf20/les_output/wrf-partmc'
# locate gridded output files in the les_output/wrf-partmc directory path
griddedoutput_files = fnmatch.filter(os.listdir(wrfpartmc_output_path), 
                                     "gridded-output_*")
if len(griddedoutput_files) != 0:
    print('Moving gridded output files')
    for file in griddedoutput_files:
        shutil.move(f"{wrfpartmc_output_path}/{file}", f"{output_path}/{file}")
