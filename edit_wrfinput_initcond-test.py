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

    src = Dataset("wrfinput_d01", 'r+')
    #dst = Dataset("wrfinput_d01_new", "w")

    # create variables if not copied over from src, modify values using checkerboard profile
    for variable_name in variables_to_modify:
        """
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
        """

        src[variable_name][0, :, :, :] = np.zeros((100, 40, 40))

        # modify ground level values
        src['so2'][0, 0, :, :] = ic_scaling*ic_array

        """
        if initcond_profile_attribs['profile_type'] == 'checkerboard_profile':
            print(f'{variable_name}: assigning checkboard')

            fx=initcond_profile_attribs['profile_fx']
            fy=initcond_profile_attribs['profile_fy']

            initcond_profile = checkerboard_profile(fx, fy, 
                                                    xgrid=xgrid, ygrid=ygrid, zgrid=zgrid, 
                                                    min_val=initcond_profile_attribs['profile_min_val'], 
                                                    max_val=initcond_profile_attribs['profile_max_val'], 
                                                    phase_shift=initcond_profile_attribs['profile_phase_shift'])
            
            # compute sum just for ground level
            initcond_profile_sum = initcond_profile[0, :, :].sum()
            if variable_name in basecase_init_totals:
                basecase_variable_total = basecase_init_totals[variable_name] # ground level sum
                if initcond_profile_sum != 0:
                    #print(basecase_variable_total, initcond_profile_sum)
                    scaling_factor = round(basecase_variable_total/initcond_profile_sum, 1) # NOTE: somehow evaluating to 1.9, debug
                else:
                    scaling_factor = 1.0
                print(f'..scaling values by {scaling_factor:3.2f}')
                # scale emission rate
                initcond_profile[0, :, :] = scaling_factor*initcond_profile[0, :, :]
                # update in json dataset
                species_initcond[variable_name]['profile_min_val'] = scaling_factor*species_initcond[variable_name]['profile_min_val']
                species_initcond[variable_name]['profile_max_val'] = scaling_factor*species_initcond[variable_name]['profile_max_val']
            else:
                print(f'..WARNING: {variable_name} initial condition not found in basecase_init_totals, IC will not be adjusted!')
    
        elif initcond_profile_attribs['profile_type'] == 'uniform':
            print(f'{variable_name}: assigning uniform profile')
            mesh = np.zeros((xgrid, ygrid, zgrid))
            mesh[:, :, :] = initcond_profile_attribs['profile_val']
            mesh = mesh.T
            initcond_profile = np.ma.array(np.array([mesh]), mask=False, dtype='float32')
        else:
            raise ValueError('invalid profile type')
        
        dst[variable_name][:] = initcond_profile
        """
        
    src.close()
    #dst.close()

    # save updates to JSON emission file
    with open('species_initcond_profile_data.json', 'w') as outfile:
        json.dump(species_initcond, outfile, indent=4)

    return

"""
def update_netcdf_names():
    os.rename("wrfinput_d01", "wrfinput_d01_unmodifed")
    os.rename("wrfinput_d01_new", "wrfinput_d01")
    return
"""

if __name__ == '__main__':
    # Access the command-line argument passed from the Bash script
    CHEM_OPT = int(sys.argv[1])
    domain_x_cells = int(sys.argv[2]) # number of grid cells in x direction
    domain_y_cells = int(sys.argv[3]) # number of grid cells in y direction
    domain_z_cells = int(sys.argv[4]) # number of grid cells in y direction

    scenario = sys.argv[5] # gas initial conditions scenario
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
    #update_netcdf_names()

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