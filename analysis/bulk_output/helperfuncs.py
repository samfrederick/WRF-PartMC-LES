from loaddatastructs import *
import numpy as np
from mcnsh import mcnormspatialhet
from nsh import normalizedspatialhet

def createEmissionsNSHDict():
    parent_dir = '/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les'
    shdir = 'spatial-het/sh-patterns'
    griddir = f'xres{domain_x_cells}yres{domain_y_cells}'

    basecase_filename = 'uniform-basecase.csv'
    basecase_array_path = os.path.join(parent_dir, shdir, griddir, basecase_filename)
    basecase_arr = np.genfromtxt(basecase_array_path, delimiter=',')

    emissions_nsh_dict = {}
    for filename in os.listdir(os.path.join(parent_dir, shdir, griddir)):
        if filename == 'fx0fy0.csv':
            continue
        scenario = filename.replace('.csv', '')
        array_path = os.path.join(parent_dir, shdir, griddir, filename)
        scenario_arr = np.genfromtxt(array_path, delimiter=',')
        scaling_factor = basecase_arr.sum() / scenario_arr.sum()
        scenario_arr = scaling_factor*scenario_arr
        #arr_nsh = mcnormspatialhet(scenario_arr, n_permutes=100000)
        arr_nsh = normalizedspatialhet(scenario_arr)
        emissions_nsh_dict[scenario] = arr_nsh
        #print(f'{scenario}, scaling factor = {scaling_factor}, NSH = {arr_nsh:4.3f}')
    emissions_nsh_dict = dict(sorted(emissions_nsh_dict.items(), key=lambda item: item[1]))

    return emissions_nsh_dict

def convertMassConctoVolumeMixingRatio(scenario_data, species_name, species_molar_weight):

    dry_air_molar_weight = 28.9644 # g/mol
    inverse_airdens = scenario_data['ALT'][:]
    MMR = inverse_airdens*scenario_data[species_name][:] # mass mixing ratio of species
    VMR = 1e6*(dry_air_molar_weight/species_molar_weight)*MMR # Volume mixing ratio in ppmv
    # need to multiply by 1e6 to convert from mol ammonium / mol air to micromol
    return VMR

def calculateNSHTimeSlice(scenario, variable):
    scenario_aerodata = aerodata_dict[scenario]['aerodata']
    
    levels = np.arange(n_levels)
    times = np.arange(n_times)
    nsh_array = np.zeros((n_times, n_levels))
    for itime in times:
        for ilevel in levels:
            level_array = scenario_aerodata[variable][itime, ilevel, :, :]
            nsh_estimate = mcnormspatialhet(level_array, n_permutes=10000)
            nsh_array[itime, ilevel] = nsh_estimate
    
    nsh_dict[scenario][variable] = nsh_array
    return nsh_array

def calculateVarZT(scenario, variable, mixingratio=True):
    levels = np.arange(n_levels)
    times = np.arange(n_times)

    var_array = np.zeros((n_times, n_levels))
    for itime in times:
        for ilevel in levels:
            if mixingratio:
                inverse_airdens = aerodata_dict[scenario]['aerodata']['ALT'][itime, ilevel, :, :]
                level_array = inverse_airdens*aerodata_dict[scenario]['aerodata'][variable][itime, ilevel, :, :]
            else:
                level_array = aerodata_dict[scenario]['aerodata'][variable][itime, ilevel, :, :]
            var_array[itime, ilevel] = level_array.mean()
    return var_array

def calculateVarPercentDiff(scenario, variable,  mixingratio=False, skip_t0=False):
    array_scenario = calculateVarZT(scenario, variable, mixingratio)
    array_basecase = calculateVarZT('basecase', variable, mixingratio)

    if skip_t0:
        array_scenario = array_scenario[1:]
        array_basecase = array_basecase[1:]

    rel_diff = 100*(array_scenario - array_basecase)/array_basecase
    return rel_diff

def calculateVarBias(scenario, variable,  mixingratio=False, skip_t0=False):
    array_scenario = calculateVarZT(scenario, variable, mixingratio)
    array_basecase = calculateVarZT('basecase', variable, mixingratio)

    if skip_t0:
        array_scenario = array_scenario[1:]
        array_basecase = array_basecase[1:]

    bias = (array_scenario - array_basecase)
    return bias

def calcTotConcZT(scenario, dist_type, min_particle_size):
    n_bins = 100
    bins = np.arange(n_bins)
    times = np.arange(n_times)

    #scenario_aerodata = aerodata_dict[scenario]['aerodata']
    scenario_distdata = aerodata_dict[scenario]['distdata']
    xgrid, ygrid, kgrid = 40, 40, 100
    var_array = np.zeros((n_times, kgrid, ygrid, xgrid))

    bin_edges = scenario_distdata['BIN_EDGES'][:].data[0]#scenario_aerodata['BIN_EDGES'][:].data[0]
    bin_centers = scenario_distdata['BIN_CENTERS'][:].data[0]#scenario_aerodata['BIN_CENTERS'][:].data[0]
    bin_width = bin_edges[1:] - bin_edges[:-1]

    for time_idx in times:
        for k in np.arange(kgrid):
            tot_conc = np.zeros((xgrid, ygrid))
            for bin_idx, width, center in zip(bins[1:], bin_width, bin_centers[1:]): 
                # just look at particles larger than 50 nanometers
                if center*1e9 < min_particle_size:
                    continue
                bin_idx += 1 # 1 indexing 
                bin_data = scenario_distdata[f'{dist_type}_a{str(bin_idx).zfill(3)}'][time_idx, k, :, :].data/1e6
                tot_conc += bin_data
            var_array[time_idx, k, :, :] = tot_conc
    var_array_mean = var_array.mean(axis=(2, 3))

    return var_array_mean

def calculateBoxplotData(variable='ccn_001'):
    boxplot_vardata = []
    for scenario_name in emissions_nsh_dict.keys():
        #scenario_nsh = emissions_nsh_dict[scenario_name]
        #scenario_nsh_values.append(scenario_nsh)
        #boxplot_positions.append(round(scenario_nsh, 2))
        #scenario_data = aerodata_dict[scenario_name]['aerodata']
        var_pdiff = calculateVarPercentDiff(scenario=scenario_name, variable=variable, mixingratio=True)
        var_pdiff = var_pdiff[18:36, :65] # previously time 24: (2hrs +)
        boxplot_vardata.append(var_pdiff.flatten())
    boxplot_data[variable] = boxplot_vardata    

emissions_nsh_dict = createEmissionsNSHDict()
scenarios = ['fx2fy2', 'fx1fy1', 'fx1fy0', 'road-double', 'point-source-10x10', 'road-8x', 'road-16x'
             ]
emissions_nsh_dict = {key: value for key, value in emissions_nsh_dict.items() if key in scenarios}