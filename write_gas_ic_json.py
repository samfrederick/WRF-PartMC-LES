import json
import sys
import os
from pathlib import Path
import pandas as pd
import json
import numpy as np

def write_gas_ic_json(job_path, gas_ic_path, scenario):
    initcond_dict = {}

    gas_species = pd.read_csv('/data/keeling/a/sf20/b/wrf_partmc/WRFV3/test/em_les/gas_params.csv', header=0)
    urban_plume_ics_emiss = pd.read_csv(gas_ic_path, header=0)

    urban_plume_gas_species = list(urban_plume_ics_emiss.Symbol.values)

    urban_plume_ics = pd.DataFrame()
    for i, species in enumerate(gas_species.Species.values):
        if species.upper() in urban_plume_gas_species:
            urban_plume_ics.loc[i, 'Species'] = species
            ic = urban_plume_ics_emiss[urban_plume_ics_emiss.Symbol==species.upper()]['Initial Mole Fraction (ppb)'].values[0]
            if np.isnan(ic):
                ic = 0
            urban_plume_ics.loc[i, 'Initial Mole Fraction (ppb)'] = ic
            
    #urban_plume_ics.to_csv('urbanplume_gas_ics.csv', index=False)

    for var, conc in zip(urban_plume_ics.Species.values, urban_plume_ics['Initial Mole Fraction (ppb)'].values):
        
        # Assign initial conditions for all anthropogenic emission variables
        # gas phase concentration values from the urban plume csv are in ppbv
        # wrfinput expects everything in ppmv so need to convert
        conc_ppmv = conc / 1000
        initcond_dict[var] = {
                        "profile_scenario": scenario,
                        "profile_min_val": 0.0,
                        "profile_max_val": conc_ppmv,
                        "profile_phase_shift": False
                        }

    with open(f'{job_path}/species_initcond_profile_data.json', 'w') as outfile:
        json.dump(initcond_dict, outfile, indent=4)


"""
def write_initcond_json(job_path, chem_opt, scenario=None):
    # initialize with all species set to near-zero values 
    init_species = read_initcond_json(chem_opt)
    init_species = set_profile_type(init_species, scenario)

    with open(f'{job_path}/species_initcond_profile_data.json', 'w') as outfile:
        json.dump(init_species, outfile, indent=4)

def read_initcond_json(chem_opt):
    # Read archival copy of profile initial condition json file
    with open(os.path.join(Path(__file__).parent, 'profile-jsons', f'chemopt{chem_opt}', 
                           'species_initcond_profile_data.json'), 'r') as infile:    
        data = json.load(infile)
    return data

def set_profile_type(species_dict, scenario):
    # set all species to emit with fx and fy frequencies
    for species in species_dict:
        species_dict[species]['profile_type'] = scenario
    return species_dict
"""

if __name__ == '__main__':
    # Access the command-line argument passed from the Bash script
    job_path = sys.argv[1]
    chem_opt = int(sys.argv[2])
    scenario = sys.argv[3]
    gas_ic_path = sys.argv[4]
    
    write_gas_ic_json(job_path, gas_ic_path, scenario)

