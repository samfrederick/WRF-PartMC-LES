"""
Generate emissions netcdf files for surface level grid cells 

Author: Sam Frederick, August 2023
"""
import sys, os
import netCDF4 as nc
import numpy as np
import pandas as pd
import shutil

from emiss_profiles import checkerboard_profile


aero_species = ['SO4','NO3','Cl','NH4','MSA','ARO1','ARO2','ALK1','OLE1','API1',
                'API2','LIM1','LIM2','CO3','Na','Ca','OIN','OC','BC','H2O']  

def load_urban_plume_aero_props():
    #TODO use os.get_cwd()
    aero_emiss = pd.read_csv('/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les/urban_plume_aero_emiss.csv')
    return aero_emiss

# ------------------------------------------------------------------------------
#
#  Routines for setting variable values
# 
# ------------------------------------------------------------------------------

def set_source_ids(n_modes, n_aero_sources):
    source_ids = []
    modes_per_source = int(n_modes/n_aero_sources)

    for i in np.arange(1, 1+n_aero_sources):
        for j in range(modes_per_source):
            source_ids.append(i)
    return np.array(source_ids)

def set_char_radii(n_modes, n_aero_sources):
    char_radii = []
    modes_per_source = int(n_modes/n_aero_sources)

    if modes_per_source == 2:
        mode_char_radii = [1.8090e-08, 3.9905e-08] # using defaults from Jeff's aero_emit_dists_001_001_001.nc file
    else:
        raise ValueError(f'Undefined mode radii for number of modes per source: {modes_per_source}')

    for i in np.arange(1, 1+n_aero_sources):
        for radius in mode_char_radii:
            char_radii.append(radius)
    return np.array(char_radii)


def set_log10_std_dev_radius(n_modes, n_aero_sources):
    log10_std_dev_radii = []
    modes_per_source = int(n_modes/n_aero_sources)

    if modes_per_source == 2:
        mode_stddev_radii = [0.20411998, 0.25527251] # using defaults from Jeff's aero_emit_dists_001_001_001.nc file
    else:
        raise ValueError(f'Undefined mode std dev for number of modes per source: {modes_per_source}')

    for i in np.arange(1, 1+n_aero_sources):
        for radius in mode_stddev_radii:
            log10_std_dev_radii.append(radius)
    return np.array(log10_std_dev_radii)

def set_num_conc(n_times, n_modes, set_zero_conc):
    # setting number concentration of each mode to zero at each time
    if set_zero_conc:
        num_conc = np.zeros((n_times, n_modes))
    else:
        print('Note: setting aerosol emissions to zero')
        num_conc = np.zeros((n_times, n_modes))
    #for i in range(n_times):
    #    num_conc[i] = np.zeros(n_modes) 
    return num_conc

def set_vol_frac(n_times, n_modes, n_aero_specs):

    vol_frac = np.zeros((n_times, n_modes, n_aero_specs))

    for itime in range(n_times):
        for imode in range(n_modes):
            for ispecies in range(n_aero_specs):
                # For now, just set the volume fraction to be equal across all species
                # Makes it easy to ensure that volume fractions add up to one
                spec_vol_frac = 1 / n_aero_specs
                vol_frac[itime, imode, ispecies] = spec_vol_frac

    return vol_frac

def set_gas_emission(n_times, n_gas_specs, set_zero_conc, emission_scaling):

    emission_rates = np.zeros((n_times, n_gas_specs))

    if set_zero_conc:
        for itime in range(n_times):
            for ispec in range(n_gas_specs):
                emission_rates[itime, ispec] = 0
    else:
        print(f'Scaling gas-phase emissions by {emission_scaling}')

        # Using emission rates from center of CARES domain, aero_emit_dists_085_080_001.nc, variable name 'gas_emission', t index = 0
        #emiss_data = pd.read_csv('/data/keeling/a/sf20/b/wrf_partmc/WRFV3/test/em_les/cares_gas_emiss_rates.csv')

        # Use emission rates from PartMC Urban-Plume scenario (rates from Riemer et al 2009)
        #TODO use os.get_cwd()
        emiss_data = pd.read_csv('/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les/urbanplume_gas_emiss_rates.csv')
        for itime in range(n_times):
            for ispec in range(n_gas_specs):
                spec_emiss_rate = emiss_data.loc[ispec, 'Emission rate (mol m^{-2} s^{-1})']
                emission_rates[itime, ispec] = emission_scaling*spec_emiss_rate

    return emission_rates


# ------------------------------------------------------------------------------
#
#  Wrapper for creating netcdf file
# 
# ------------------------------------------------------------------------------

def create_emit_ncfile(set_zero_conc=False, emission_scaling=1.0):

    if set_zero_conc:
        file_prefix = 'nonzero'
    else:
        file_prefix = 'zero'
        
    filename = f'aeroemit_{file_prefix}.nc'

    print(f'Scaling aerosol emissions by {emission_scaling}')

    ncfile = nc.Dataset(filename, 'w')

    URBAN_PLUME_EMISS = True

    # Dimension values
    n_times = 24
    n_modes = 3 # 1 aerosol mode per urban plume source
    n_aero_specs = 20
    n_aero_sources = 3 #3 aerosol sources for urban plume (cooking, diesel, gas)
    n_gas_specs = 77

    one_hour = 3600

    # Create dimensions
    times_dim = ncfile.createDimension('n_times', n_times)
    modes_dim = ncfile.createDimension('n_modes', n_modes)
    aero_specs_dim = ncfile.createDimension('n_aero_specs', n_aero_specs)
    aero_sources_dim = ncfile.createDimension('n_aero_sources', n_aero_sources)
    gas_specs_dim = ncfile.createDimension('n_gas_specs', n_gas_specs)

    # Create variables, add descriptions
    aero_emission_rate_scale = ncfile.createVariable('aero_emission_rate_scale', 'f8', ('n_times'))
    aero_emission_rate_scale.unit = '(1)'
    aero_emission_rate_scale.long_name = "Aerosol emission rate"
    aero_emission_rate_scale.description = "Aerosol emission rate scales at set-points"
    
    aero_emission_time = ncfile.createVariable('aero_emission_time', 'f8', ('n_times'))
    aero_emission_time.unit = 's'
    aero_emission_time.long_name = "Aerosol emission time"
    aero_emission_time.description = "Aerosol emission set-points times (s)."
    
    char_radius = ncfile.createVariable('char_radius', 'f8', ('n_modes'))
    char_radius.unit = 'm'
    char_radius.long_name = "characteristic_radius"
    char_radius.standard_name = "characteristic_radius"
    char_radius.description = "Characteristic radius, with meaning dependent on mode type"
    
    log10_std_dev_radius = ncfile.createVariable('log10_std_dev_radius', 'f8', ('n_modes'))
    log10_std_dev_radius.unit = 'm'
    log10_std_dev_radius.long_name = "log10_std_dev_radius"
    log10_std_dev_radius.standard_name = "log10_std_dev_radius"
    log10_std_dev_radius.description = "Log base 10 of geometric standard deviation of radius, (m)."
    
    source_id = ncfile.createVariable('source_id', 'i', ('n_modes'))
    source_id.unit = '(1)'
    source_id.long_name = "Source number."
    source_id.standard_name = "Source number"
    source_id.description = "Source number for each emission mode."
    
    num_conc = ncfile.createVariable('num_conc', 'f8', ('n_times', 'n_modes'))
    num_conc.unit = '# m^{-3}'
    num_conc.long_name = "total number concentration"
    num_conc.standard_name = "total number concentration"
    num_conc.description = "Total number concentration of mode (#/m^3)."
    
    vol_frac = ncfile.createVariable('vol_frac', 'f8', ('n_times', 'n_modes', 'n_aero_specs'))
    vol_frac.unit = '(1)'
    vol_frac.long_name = "species fractions"
    vol_frac.standard_name = "species_fractions"
    vol_frac.description = "Species fractions by volume [length \\c aero_data%%n_spec]."

    source_weight_class = ncfile.createVariable('source_weight_class', 'i', ('n_modes'))
    source_weight_class.unit = '(1)'
    source_weight_class.long_name = "Aerosol weight class"
    source_weight_class.description = "Weight class ID for each aerosol mode."

    aero_source = ncfile.createVariable('aero_source', 'i', ('n_aero_sources'))
    aero_source.names = 'Meat-cooking,Diesel-vehicles,Gasoline-vehicles' # note the number of comma-delimited source names should match the n_aero_sources dimension
    aero_source.description = "dummy dimension variable (no useful value) - read source names as comma-separated values from the \'names\' attribute"

    # Set aerosol emission values
    if URBAN_PLUME_EMISS:
        aero_props = load_urban_plume_aero_props()
        aero_emission_rate_scale[:] = 1 # set all rates to 1
        aero_emission_time[:] = np.arange(0, n_times*one_hour, one_hour) # set times to each hour during a 24hr period

        # set the characteristic radius of each mode (geometric mean radius)
        char_radius[:] = np.array(0.5*aero_props.geom_mean_diam)

        # set the log10 of std dev of radius (equivalent to log10 of geom std dev of diameter)
        log10_std_dev_radius[:] = np.array(aero_props.log10_geom_std_dev)

        # set source ids
        source_id[:] = set_source_ids(n_modes, n_aero_sources)

        # set the number concentration
        num_conc[:] = np.zeros((n_times, n_modes))
        if not set_zero_conc:
            #num_conc[:] = emission_scaling*aero_props.num_conc

            # Set aerosol num conc to zero for the first hour of the simulation
            # NOTE that the timing in the aero_emit_dist is relative to the duration 
            # of the simulation rather than UTC time
            num_conc[0, :] = 0
            num_conc[1:, :] = emission_scaling*aero_props.num_conc

        # set the volume fraction of species
        OC_idx = 17
        BC_idx = 18
        cooking_mode_idx = 0
        diesel_mode_idx = 1
        gas_mode_idx = 2
        vol_frac[:] = 0 # Bugfix: set all values zero intiiallly (otherwise values are nans!)
        vol_frac[:, cooking_mode_idx, OC_idx] = 1
        vol_frac[:, diesel_mode_idx, OC_idx] = .3
        vol_frac[:, diesel_mode_idx, BC_idx] = .7
        vol_frac[:, gas_mode_idx, OC_idx] = .8
        vol_frac[:, gas_mode_idx, BC_idx] = .2

        source_weight_class[:] = 1 # just set all wieght classes to 1 for now
        aero_source[:] = np.arange(1, n_aero_sources+1, 1)
    else:
        aero_emission_rate_scale[:] = 1 # set all rates to 1
        aero_emission_time[:] = np.arange(0, n_times*one_hour, one_hour) # set times to each hour during a 24hr period
        char_radius[:] = set_char_radii(n_modes, n_aero_sources)
        log10_std_dev_radius[:] = set_log10_std_dev_radius(n_modes, n_aero_sources)
        source_id[:] = set_source_ids(n_modes, n_aero_sources)
        num_conc[:] = set_num_conc(n_times, n_modes, set_zero_conc)
        vol_frac[:] = set_vol_frac(n_times, n_modes, n_aero_specs)
        source_weight_class[:] = 1 # just set all wieght classes to 1 for now
        aero_source[:] = np.arange(1, n_aero_sources+1, 1)

    # NOTE: Jeff has these values set to 5e-2 for most gas species for some reason?
    gas_emission_rate_scale = ncfile.createVariable('gas_emission_rate_scale', 'f8', ('n_times'))
    gas_emission_rate_scale.unit = '(1)'
    gas_emission_rate_scale.long_name = "Gas emission rate scale factor"
    gas_emission_rate_scale.description = "Gas emission rate scales at set-points"
    gas_emission_rate_scale[:] = 1 # set all rates to 1

    gas_emission_time = ncfile.createVariable('gas_emission_time', 'f8', ('n_times'))
    gas_emission_time.unit = 's'
    gas_emission_time.long_name = "Gas emission time"
    gas_emission_time.description = "Gas emission set-points times (s)."
    gas_emission_time[:] = np.arange(0, n_times*one_hour, one_hour) # set times to each hour during a 24hr period

    gas_emission = ncfile.createVariable('gas_emission', 'f8', ('n_times', 'n_gas_specs'))
    gas_emission.unit = 'mol m^{-2} s^{-1}'
    gas_emission.long_name = "gas emissions"
    gas_emission.standard_name = "gas emissions"
    gas_emission.description = "gas phase emission rates."
    gas_emission[:] = set_gas_emission(n_times, n_gas_specs, set_zero_conc, emission_scaling)

    # Set the gas emissions rates to zero for first hour of the simulation
    if URBAN_PLUME_EMISS:
        gas_emission[0, :] = 0
    
    ncfile.close()

    return filename

if __name__ == '__main__':

    scenario = sys.argv[1] # emissions scenario

    domain_x_cells = int(sys.argv[2]) # number of grid cells in x direction
    domain_y_cells = int(sys.argv[3]) # number of grid cells in y direction
    domain_z_cells = int(sys.argv[4]) # number of grid cells in y direction

    #TODO use os.get_cwd()
    cwd = '/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les'
    aero_emit_dist_path = '/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les/aero_emit_dists'

    print(f'Setting aerosol and gas emissions with scenario: {scenario}')

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

    filename_nonzero_emiss = create_emit_ncfile(set_zero_conc=False, 
                                                emission_scaling=scaling_factor)
    filename_zero_emiss = create_emit_ncfile(set_zero_conc=True, 
                                             emission_scaling=scaling_factor)

    shutil.move(os.path.join(cwd, filename_nonzero_emiss), 
                    os.path.join(aero_emit_dist_path,filename_nonzero_emiss))
    shutil.move(os.path.join(cwd, filename_zero_emiss), 
                    os.path.join(aero_emit_dist_path,filename_zero_emiss))

    # These are dummy, binary values to determine whether emissions are present 
    # in an area or not (not emission rates)
    #max_val = 1
    #min_val = 0
    #checkerboard_z_cells = 1 # number of grid cells in z direction (arbitrary if just concerned with ground pattern)
    #checkerboard = checkerboard_profile(fx, fy, domain_x_cells+1, domain_y_cells+1, 
    #                                    checkerboard_z_cells, max_val, min_val, phase_shift=False)[0][0]

    for i in np.arange(1, domain_x_cells+1):
        for j in np.arange(1, domain_y_cells+1):

            # for some reason I need to reverse the indexing to match the domain
            # initialization for gas-phase species
            grid_value = scenario_arr[j-1, i-1] 
            str_i = str(i).zfill(3)
            str_j = str(j).zfill(3)
            k = '001'
            destination_filename = f'aero_emit_dists_{str_i}_{str_j}_{k}.nc'

            if grid_value == 1.0:
                #print('non-zero emissions:', i, j)
                shutil.copyfile(os.path.join(aero_emit_dist_path, filename_nonzero_emiss), os.path.join(aero_emit_dist_path,destination_filename))
            else:
                #print('no emissions:', i, j)
                shutil.copyfile(os.path.join(aero_emit_dist_path, filename_zero_emiss), os.path.join(aero_emit_dist_path, destination_filename))
            
            #print(destination_filename)