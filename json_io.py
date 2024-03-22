import json
import sys
import os
from pathlib import Path

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

if __name__ == '__main__':
    # Access the command-line argument passed from the Bash script
    job_path = sys.argv[1]
    chem_opt = int(sys.argv[2])
    scenario = sys.argv[3]
    
    write_initcond_json(job_path, chem_opt, scenario=scenario)#overlap_precursors=overlap_precursors)

