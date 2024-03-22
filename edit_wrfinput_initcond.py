"""
"""
import os#, sys
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
#from emiss_profiles import gauss3d_profile, checkerboard_profile
import json
import sys

def create_modified_netcdf(job_path, ic_array=None, ic_scaling=1.0):
    # Code via 
    # https://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one
    # User Rich Signell

    # Modified by Sam Frederick, October 2022

    jsonfile = open(f'{job_path}/species_initcond_profile_data.json')
    species_initcond = json.load(jsonfile)
    variables_to_modify = [species for species in species_initcond.keys()]

    src = Dataset(f"{job_path}/wrfinput_d01")
    dst = Dataset(f"{job_path}/wrfinput_d01_new", "w")

    # Mesh grid dimensions
    xgrid = src.dimensions['west_east'].size
    ygrid = src.dimensions['south_north'].size
    zgrid = src.dimensions['bottom_top'].size

    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for variable_name, dimension in src.dimensions.items():
        dst.createDimension(
            variable_name, (len(dimension) if not dimension.isunlimited() else None))

    # copy over all variables from src to dst
    for variable_name, variable_attribs in src.variables.items():
        #print(f'..copying over data for {variable_name}')
        x = dst.createVariable(variable_name, variable_attribs.datatype, variable_attribs.dimensions)
        dst[variable_name][:] = src[variable_name][:]
        # copy variable attributes all at once via dictionary
        # NOTE: this assigns "MemoryOrder", "description", "units", "stagger", and "fieldtype"
        dst[variable_name].setncatts(src[variable_name].__dict__)

    # create variables if not copied over from src, modify values using checkerboard profile
    for variable_name in variables_to_modify:
        if variable_name.lower() not in [var.lower() for var in src.variables.keys()]:

            print(f'..{variable_name}: not present in source wrfinput file, creating variable using so2 attributes')
            reference_variable = src.variables['so2']
            x = dst.createVariable(variable_name, reference_variable.datatype, reference_variable.dimensions)
            #  NOTE: this assigns "MemoryOrder", "description", "units", "stagger", and "fieldtype"
            dst[variable_name].setncatts(src[reference_variable.name].__dict__)

        cell_conc = species_initcond[variable_name]['profile_max_val']
        dst[variable_name][0, :, :, :] = np.zeros((zgrid, ygrid, xgrid))
        # modify ground level values
        #dst[variable_name][0, 0, :, :] = ic_scaling*cell_conc*ic_array
        print(f'{variable_name} concentration:', ic_scaling*cell_conc*ic_array.max())
        
        # EDIT 2/16/24 - set emissions to uniform in domain
        dst[variable_name][0, :, :, :] = ic_scaling*cell_conc*ic_array

    print('\n')
    src.close()
    dst.close()

    # save updates to JSON emission file
    with open(f'{job_path}/species_initcond_profile_data.json', 'w') as outfile:
        json.dump(species_initcond, outfile, indent=4)

    return

def update_netcdf_names(job_path):
    os.rename(f'{job_path}/wrfinput_d01', f'{job_path}/wrfinput_d01_unmodifed')
    os.rename(f'{job_path}/wrfinput_d01_new', f'{job_path}/wrfinput_d01')
    return

if __name__ == '__main__':
    # Access the command-line argument passed from the Bash script
    job_path = sys.argv[1]
    CHEM_OPT = int(sys.argv[2])
    domain_x_cells = int(sys.argv[3]) # number of grid cells in x direction
    domain_y_cells = int(sys.argv[4]) # number of grid cells in y direction
    domain_z_cells = int(sys.argv[5]) # number of grid cells in y direction

    #scenario = sys.argv[5] # gas initial conditions scenario
    # EDIT 2/16/24 - set initial condition to be uniform
    scenario = 'uniform-basecase'
    print('\n')
    print(f'Setting gas-phase initial conditions with scenario: {scenario}')

    simdir = Path(__file__).parent
    shdir = 'spatial-het/sh-patterns'
    griddir = f'xres{domain_x_cells}yres{domain_y_cells}'
    filename = f'{scenario}.csv'
    array_path = os.path.join(simdir, shdir, griddir, filename)
    scenario_arr = np.genfromtxt(array_path, delimiter=',')

    basecase_filename = 'uniform-basecase.csv'
    basecase_array_path = os.path.join(simdir, shdir, griddir, basecase_filename)
    basecase_arr = np.genfromtxt(basecase_array_path, delimiter=',')

    scaling_factor = basecase_arr.sum() / scenario_arr.sum()

    create_modified_netcdf(job_path, ic_array=scenario_arr, ic_scaling=scaling_factor)
    update_netcdf_names(job_path)
    