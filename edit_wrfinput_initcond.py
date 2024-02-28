"""
"""
import os#, sys
import numpy as np
from netCDF4 import Dataset
#from emiss_profiles import gauss3d_profile, checkerboard_profile
import json
import sys

def create_modified_netcdf(ic_array=None, ic_scaling=1.0):
    # Code via 
    # https://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one
    # User Rich Signell

    # Modified by Sam Frederick, October 2022

    jsonfile = open('species_initcond_profile_data.json')
    species_initcond = json.load(jsonfile)
    variables_to_modify = [species for species in species_initcond.keys()]

    src = Dataset("wrfinput_d01")
    dst = Dataset("wrfinput_d01_new", "w")

    # Mesh grid dimensions
    xgrid = src.dimensions['west_east'].size
    ygrid = src.dimensions['south_north'].size
    zgrid = src.dimensions['bottom_top'].size
    """
    
    # species initial conidtion concentrations summed over the entire domain for the 
    # basecase where species initialized everywhere with concentration indicated
    # in species_initcond_profile_data.json
    basecase_init_totals = {}
    for species in species_initcond:
        if species_initcond[species]['profile_type'] == 'checkerboard_profile':
            per_cell_conc = species_initcond[species]['profile_max_val']
            sum_ground_conc = per_cell_conc*(xgrid)*(ygrid)  
            basecase_init_totals[species] = sum_ground_conc
    """

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
            
            #if variable_name in aero_species['aerosol_species[ug/kg-dryair]'].values:
            #    print(f'{variable_name}: aerosol species, setting units to "ug/kg-dryair"')
            #    dst[variable_name].units = "ug/kg-dryair"
        
        #initcond_profile_attribs = species_initcond[variable_name]

        cell_conc = species_initcond[variable_name]['profile_max_val']
        dst[variable_name][0, :, :, :] = np.zeros((zgrid, ygrid, xgrid))
         # modify ground level values
        print(f'{variable_name} ground concentration:', ic_scaling*cell_conc*ic_array.max())
        #dst[variable_name][0, 0, :, :] = ic_scaling*cell_conc*ic_array
        
        # EDIT 2/16/24 - set emissions to uniform in domain
        dst[variable_name][0, :, :, :] = ic_scaling*cell_conc*ic_array

    print('\n')
    src.close()
    dst.close()

    # save updates to JSON emission file
    with open('species_initcond_profile_data.json', 'w') as outfile:
        json.dump(species_initcond, outfile, indent=4)

    return

def update_netcdf_names():
    os.rename("wrfinput_d01", "wrfinput_d01_unmodifed")
    os.rename("wrfinput_d01_new", "wrfinput_d01")
    return

if __name__ == '__main__':
    # Access the command-line argument passed from the Bash script
    CHEM_OPT = int(sys.argv[1])
    domain_x_cells = int(sys.argv[2]) # number of grid cells in x direction
    domain_y_cells = int(sys.argv[3]) # number of grid cells in y direction
    domain_z_cells = int(sys.argv[4]) # number of grid cells in y direction

    #scenario = sys.argv[5] # gas initial conditions scenario
    # EDIT 2/16/24 - set initial condition to be uniform
    scenario = 'uniform-basecase'
    print('\n')
    print(f'Setting gas-phase initial conditions with scenario: {scenario}')

    cwd = os.getcwd()
    shdir = 'spatial-het/sh-patterns'
    griddir = f'xres{domain_x_cells}yres{domain_y_cells}'
    filename = f'{scenario}.csv'
    array_path = os.path.join(cwd, shdir, griddir, filename)
    scenario_arr = np.genfromtxt(array_path, delimiter=',')

    basecase_filename = 'uniform-basecase.csv'
    basecase_array_path = os.path.join(cwd, shdir, griddir, basecase_filename)
    basecase_arr = np.genfromtxt(basecase_array_path, delimiter=',')

    scaling_factor = basecase_arr.sum() / scenario_arr.sum()

    create_modified_netcdf(ic_array=scenario_arr, ic_scaling=scaling_factor)
    update_netcdf_names()

    # a list of the species in WRF PartMC if needed...
    species = ['h2so4', 'hno3', 'hcl', 'nh3', 'no', 'no2', 'no3', 'n2o5', 
               'hono', 'hno4', 'o3', 'o1d', 'O3P', 'oh', 'ho2', 'h2o2', 'co',
               'so2', 'ch4', 'c2h6', 'ch3o2', 'ethp', 'hcho', 'ch3oh', 'ANOL',
               'ch3ooh', 'ETHOOH', 'ald2', 'hcooh', 'RCOOH', 'c2o3', 'pan', 
               'aro1', 'aro2', 'alk1', 'ole1', 'api1', 'api2', 'lim1', 'lim2',
               'par', 'AONE', 'mgly', 'eth', 'OLET', 'OLEI', 'tol', 'xyl', 
               'cres', 'to2', 'cro', 'open', 'onit', 'rooh', 'ro2', 'ano2', 
               'nap', 'xo2', 'xpar', 'isop', 'isoprd', 'isopp', 'isopn', 'isopo2', 
               'api', 'lim', 'dms', 'msa', 'dmso', 'dmso2', 'ch3so2h', 'ch3sch2oo', 
               'ch3so2', 'ch3so3', 'ch3so2oo', 'ch3so2ch2oo', 'SULFHOX']