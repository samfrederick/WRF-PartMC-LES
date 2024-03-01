import pandas as pd
import json
import numpy as np


initcond_dict = {}

gas_species = pd.read_csv('/data/keeling/a/sf20/b/wrf_partmc/WRFV3/test/em_les/gas_params.csv', header=0)
urban_plume_ics_emiss = pd.read_csv('/data/keeling/a/sf20/b/wrf_partmc/WRFV3/test/em_les/urban-plume-gas-emissions-ics.csv', header=0)

urban_plume_gas_species = list(urban_plume_ics_emiss.Symbol.values)

urban_plume_ics = pd.DataFrame()
for i, species in enumerate(gas_species.Species.values):
    if species.upper() in urban_plume_gas_species:
        urban_plume_ics.loc[i, 'Species'] = species
        ic = urban_plume_ics_emiss[urban_plume_ics_emiss.Symbol==species.upper()]['Initial Mole Fraction (ppb)'].values[0]
        if np.isnan(ic):
              ic = 0
        urban_plume_ics.loc[i, 'Initial Mole Fraction (ppb)'] = ic

        
urban_plume_ics.to_csv('urbanplume_gas_ics.csv', index=False)

for var, conc in zip(urban_plume_ics.Species.values, urban_plume_ics['Initial Mole Fraction (ppb)'].values):
      
       # Assign initial conditions for all anthropogenic emission variables
       # gas phase concentration values are in ppbv
       # wrfinput expects everything in ppmv so need to convert
       conc_ppmv = conc / 1000
       initcond_dict[var] = {
                     #"profile_type": "checkerboard_profile",
                     #"profile_fx": 0.5,
                     #"profile_fy": 0.5,
                     "profile_scenario": "uniform-basecase",
                     "profile_min_val": 0.0,
                     "profile_max_val": conc_ppmv,
                     "profile_phase_shift": False
                     }

with open(f'species_initcond_profile_data.json', 'w') as outfile:
        json.dump(initcond_dict, outfile, indent=4)
