import os
import numpy as np
import pandas as pd

USE_EXACT_SH = True

if USE_EXACT_SH:
    from nsh import normalizedspatialhet
    file_suffix = 'exact'
else:
    from mcnsh import mcnormspatialhet
    file_suffix = 'mc-estimate'

domain_x_cells = 100
domain_y_cells = 100

scenario_df = pd.DataFrame(columns=['scenario', 'NSH', 'scaling-factor'])

cwd = os.getcwd()
shdir = 'sh-patterns'
griddir = f'xres{domain_x_cells}yres{domain_y_cells}'

basecase_filename = 'uniform-basecase.csv'
basecase_array_path = os.path.join(cwd, shdir, griddir, basecase_filename)
basecase_arr = np.genfromtxt(basecase_array_path, delimiter=',')

nsh_dict = {}
i = 0
for filename in os.listdir(os.path.join(cwd, shdir, griddir)):
    if filename == 'fx0fy0.csv':
        continue
    
    scenario = filename.replace('.csv', '')
    print(f'{scenario}')
    array_path = os.path.join(cwd, shdir, griddir, filename)
    scenario_arr = np.genfromtxt(array_path, delimiter=',')
    scaling_factor = basecase_arr.sum() / scenario_arr.sum()
    scenario_arr = scaling_factor*scenario_arr

    if USE_EXACT_SH:
        arr_nsh = normalizedspatialhet(scenario_arr)
    else:
        arr_nsh = mcnormspatialhet(scenario_arr, n_permutes=100000)
    
    nsh_dict[scenario] = arr_nsh
    print(f'..scaling factor = {scaling_factor}, NSH = {arr_nsh:4.3f}')

    scenario_df.loc[i, 'scenario'] = scenario
    scenario_df.loc[i, 'NSH'] = np.round(arr_nsh, 4)
    scenario_df.loc[i, 'scaling-factor'] = np.round(scaling_factor, 2)
    i += 1

scenario_df = scenario_df.sort_values(by='NSH')
scenario_df.to_csv(f'sh_patterns_xres{domain_x_cells}_yres{domain_y_cells}_{file_suffix}.csv', index=False)
#sorted_nsh_dict = dict(sorted(nsh_dict.items(), key=lambda item: item[1]))
