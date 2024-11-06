"""
Create initial condition netcdf files for wrf-partmc.

Presently this simply copies one of the data files in the center of the CARES 
domain and uses this across the entire domain.
"""
import shutil
from pathlib import Path
import numpy as np
import os, sys
import json
import netCDF4 as nc
from emiss_profiles import checkerboard_profile

"""
urban_plume_aero_ics = {'aitken_mode': {'number_conc [m^-3]': 3.2e9, # 3200 per cubic cm
                                        'geom_mean_diam [m]': 0.02e-6,
                                        'geom_std_dev': 1.45,
                                        'mass_frac': {'SO4': 0.3636,'NO3': 0,
                                                      'Cl': 0,'NH4': 0.1364,
                                                      'MSA': 0,'ARO1': 0,
                                                      'ARO2': 0, 'ALK1': 0,
                                                      'OLE1': 0,'API1': 0,
                                                      'API2': 0,'LIM1': 0,
                                                      'LIM2': 0,'CO3': 0,
                                                      'Na': 0,'Ca': 0,
                                                      'OIN': 0,'OC': 0.5,
                                                      'BC': 0,'H2O': 0}
                                         },
                        'accum_mode': {'number_conc [m^-3]': 2.9e9,
                                        'geom_mean_diam [m]': 0.116e-6,
                                        'geom_std_dev': 1.65,
                                        'mass_frac': {'SO4': 0.3636,'NO3': 0,
                                                      'Cl': 0,'NH4': 0.1364,
                                                      'MSA': 0,'ARO1': 0,
                                                      'ARO2': 0, 'ALK1': 0,
                                                      'OLE1': 0,'API1': 0,
                                                      'API2': 0,'LIM1': 0,
                                                      'LIM2': 0,'CO3': 0,
                                                      'Na': 0,'Ca': 0,
                                                      'OIN': 0,'OC': 0.5,
                                                      'BC': 0,'H2O': 0}
                                         },
                        }
"""
def mass_frac_to_volume_frac(mass_frac, species_density):
    return (mass_frac/species_density)/((mass_frac/species_density).sum())

def createICReference(job_path, aero_ics_path, set_zero_conc=False, domain_z_cells=None, ic_scaling=1.0):

    if set_zero_conc:
        file_prefix = 'nonzero'
    else:
        file_prefix = 'zero'
    
    with open(aero_ics_path, 'r') as infile:    
        urban_plume_aero_ics = json.load(infile)

    with open(os.path.join(os.path.dirname(aero_ics_path), 'species-densities.json'), 'r') as infile:
        species_density = json.load(infile)
        species_density = np.array(list(species_density.values()))
    
    file_path = f'{job_path}/ics/ic_{file_prefix}_reference.nc'
    ncfile = nc.Dataset(file_path, 'w')

    if domain_z_cells is None:
        raise AttributeError('Number of vertical grid cells must be specified')
    
    print(f'Scaling aerosol initial conditions by {ic_scaling}')

    n_levels = domain_z_cells
    n_aero_modes = 2
    n_aero_specs = 20

    modes_dim = ncfile.createDimension('n_aero_modes', n_aero_modes)
    aero_specs_dim = ncfile.createDimension('n_aero_specs', n_aero_specs)
    z_dim = ncfile.createDimension('nz', n_levels)

    # TODO: not sure what this is for...Jeff has this set to 1 for all modes, levels
    mode_type = ncfile.createVariable('mode_type', 'f8', ('n_aero_modes', 'nz'))
    mode_type[:] = 1

    # float64 char_radius(n_aero_modes)
    char_radius = ncfile.createVariable('char_radius', 'f8', ('n_aero_modes', 'nz'))
    char_radius.unit = 'm'
    char_radius.long_name = "characteristic_radius"
    char_radius.standard_name = "characteristic_radius"
    char_radius.description = "Characteristic radius, with meaning dependent on mode type"
    char_radii = []
    for mode, attributes in urban_plume_aero_ics.items():
        geom_mean_diam = attributes['geom_mean_diam [m]']
        geom_mean_radius = 0.5*geom_mean_diam
        char_radii.append(geom_mean_radius)
    char_radii = np.repeat(np.array(char_radii).reshape(n_aero_modes,1), repeats=n_levels, axis=1)
    char_radius[:, :] = char_radii

    # float64 log10_std_dev_radius(n_aero_modes)
    log10_std_dev_radius = ncfile.createVariable('log10_std_dev_radius', 'f8', ('n_aero_modes', 'nz'))
    log10_std_dev_radius.unit = 'm'
    log10_std_dev_radius.long_name = "log10_std_dev_radius"
    log10_std_dev_radius.standard_name = "log10_std_dev_radius"
    log10_std_dev_radius.description = "Log base 10 of geometric standard deviation of radius, (m)."
    log10_std_dev_radii = []
    for mode, attributes in urban_plume_aero_ics.items():
        std_dev = attributes['geom_std_dev']
        log10_std_dev = np.log10(std_dev)
        log10_std_dev_radii.append(log10_std_dev)
    log10_std_dev_radii = np.repeat(np.array(log10_std_dev_radii).reshape(n_aero_modes,1), repeats=n_levels, axis=1)
    log10_std_dev_radius[:] = log10_std_dev_radii

    # float64 num_conc(n_aero_modes)
    num_conc = ncfile.createVariable('num_conc', 'f8', ('n_aero_modes', 'nz'))
    num_conc.unit = '# m^{-3}'
    num_conc.long_name = "total number concentration"
    num_conc.standard_name = "total number concentration"
    num_conc.description = "Total number concentration of mode (#/m^3)."
    num_concentrations = []
    if set_zero_conc:
        num_conc[:] = 0
    else:
        for mode, attributes in urban_plume_aero_ics.items():
            num_concentration = attributes['number_conc [m^-3]']
            num_concentrations.append(num_concentration)
        num_concentrations = np.repeat(np.array(num_concentrations).reshape(n_aero_modes,1), repeats=n_levels, axis=1)
        num_conc[:] = ic_scaling*num_concentrations

    # float64 vol_frac(n_aero_modes, n_aero_specs)
    # float64 vol_frac(n_modes, n_aero_specs)
    vol_frac = ncfile.createVariable('vol_frac', 'f8', ('n_aero_specs', 'n_aero_modes', 'nz'))
    vol_frac.unit = '(1)'
    vol_frac.long_name = "species fractions"
    vol_frac.standard_name = "species_fractions"
    vol_frac.description = "Species fractions by volume [length \c aero_data%%n_spec]."
    mass_frac_array = np.zeros((n_aero_specs, n_aero_modes))
    for i, (mode, attributes) in enumerate(urban_plume_aero_ics.items()):
        mass_fracs = attributes['mass_frac'] # could convert to list
        mass_frac_list = list(mass_fracs.values())
        mass_frac_array[:, i] = np.array(mass_frac_list)
    vol_frac_array = np.zeros((n_aero_specs, n_aero_modes, n_levels))
    """
    for i in range(n_levels):
        # NOTE: IMPORTANT! This is INCORRECT - need to convert species mass fraction to volume fraction here
        # Convert as vol_frac_array = (mass_fraction_array/species_density)/((mass_fraction_array/species_density).sum())
        # where vol_frac_array is a (20,) array of volume fracs, species density is a (20,) array of species densities 
        # and mass_fraction_array is a (20,) array of species mass fractions
        #vol_frac_array = (mass_frac_array/species_density)/((mass_frac_array/species_density).sum())
        vol_frac_array[:, :, i] = mass_frac_array
    """
    for i in range(n_levels):
        for j in range(n_aero_modes):
            vol_frac_array[:,j,i] = mass_frac_to_volume_frac(mass_frac_array[:,j], species_density)
    vol_frac[:] = vol_frac_array

    vol_frac_std = ncfile.createVariable('vol_frac_std', 'f8', ('n_aero_specs', 'n_aero_modes', 'nz'))
    vol_frac_std.unit = '(1)'
    vol_frac_std.long_name = "species fractions"
    vol_frac_std.standard_name = "species_fractions"
    vol_frac_std.description = "Species fractions by volume [length \c aero_data%%n_spec]."
    #mass_frac_array = np.zeros((n_aero_modes, n_aero_specs))
    #for i, (mode, attributes) in enumerate(urban_plume_aero_ics.items()):
    #    mass_fracs = attributes['mass_frac'] # could convert to list
    #    mass_frac_list = list(mass_fracs.values())
    #    mass_frac_array[i, :] = np.array(mass_frac_list)
    vol_frac_std[:] = 0 # Jeff seems to have all of these set to zero

    # int32 source(n_aero_modes) 
    # NOTE: this is source_id in the emission files
    source = ncfile.createVariable('source', 'i', ('n_aero_modes', 'nz'))
    source.unit = '(1)'
    source.long_name = "Source number."
    source.standard_name = "Source number"
    source.description = "Source number."
    sources = np.repeat(np.arange(1, 3, 1).reshape(n_aero_modes,1), repeats=n_levels, axis=1)
    source[:] = sources

    ncfile.close()

    return file_path

if __name__ == '__main__':
    

    # SF edit 11/9/23 - change aerosol ICs to uniform over the entire domain
    job_path = sys.argv[1]
    scenario = sys.argv[2] # initial condition scenario
    #scenario = 'uniform-basecase'

    domain_x_cells = int(sys.argv[3]) # number of grid cells in x direction
    domain_y_cells = int(sys.argv[4]) # number of grid cells in y direction
    domain_z_cells = int(sys.argv[5]) # number of grid cells in y direction

    aero_ics_path = sys.argv[6]

    print(f'Setting aerosol ICs with scenario: {scenario}')

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

    filepath_nonzero_ic = createICReference(job_path,
                                            aero_ics_path,
                                            set_zero_conc=False,
                                            domain_z_cells=domain_z_cells, 
                                            ic_scaling=scaling_factor)
    filepath_zero_ic = createICReference(job_path,
                                         aero_ics_path,
                                         set_zero_conc=True, 
                                         domain_z_cells=domain_z_cells, 
                                         ic_scaling=scaling_factor)

    # Simply copy the data from the first netcdf file into netcdf files for all gridcells
    for i in np.arange(1, domain_x_cells+1):
        for j in np.arange(1, domain_y_cells+1):

            # for some reason I need to reverse the indexing to match the domain
            # initialization for gas-phase species
            grid_value = scenario_arr[j-1, i-1] 
            str_i = str(i).zfill(3)
            str_j = str(j).zfill(3)

            destination_filename = f'ics_{str_i}_{str_j}.nc'
            destination_file_path = os.path.join(f'{job_path}/ics', destination_filename)

            if grid_value == 1.0:
                #print('non-zero emissions:', i, j)
                shutil.copyfile(filepath_nonzero_ic, destination_file_path)
            else:
                #print('no emissions:', i, j)
                shutil.copyfile(filepath_zero_ic, destination_file_path)
