import os
import netCDF4 as nc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors

class DataStruct:
    wrf_vars = ['T', 'P', 'ALT', 'PB', 'DNW', 'DN', 'Z', 'Z_AT_W', 'MAPFAC_M', 'MAPFAC_U', 'MAPFAC_V', 'MAPFAC_MX', 
            'MAPFAC_MY', 'MAPFAC_UX', 'MAPFAC_UY', 'MAPFAC_VX', 'MF_VX_INV', 'MAPFAC_VY','DENSITY_DRY_AIR', 
            'TEMPERATURE', 'REL_HUMID', 'XLAT', 'XLONG', 'ZNU', 'ZNW', 'U', 'MU_U', 'V', 'MU_V', 'W', 'WW', 'RW']
            
    aero_vars = ['TAUAER1', 'TAUAER2', 'TAUAER3', 'TAUAER4', 'GAER1', 'GAER2', 'GAER3', 'GAER4', 'WAER1', 'WAER2', 
                'WAER3', 'WAER4', 'NUM_CONC_a01', 'NUM_CONC_a02', 'NUM_CONC_a03', 'NUM_CONC_a04', 'NUM_CONC_a05',
                'NUM_CONC_a06', 'NUM_CONC_a07', 'NUM_CONC_a08', 'NUM_CONC_a09', 'NUM_CONC_a10', 'NUM_CONC_a11', 
                'NUM_CONC_a12', 'NUM_CONC_a13', 'NUM_CONC_a14', 'NUM_CONC_a15', 'NUM_CONC_a16', 'NUM_CONC_a17', 
                'NUM_CONC_a18', 'NUM_CONC_a19', 'NUM_CONC_a20', 'NUM_CONC_a21', 'NUM_CONC_a22', 'NUM_CONC_a23', 
                'NUM_CONC_a24', 'NUM_CONC_a25', 'NUM_CONC_a26', 'NUM_CONC_a27', 'NUM_CONC_a28', 'NUM_CONC_a29', 
                'NUM_CONC_a30', 'NUM_CONC_a31', 'NUM_CONC_a32', 'NUM_CONC_a33', 'NUM_CONC_a34', 'NUM_CONC_a35', 
                'NUM_CONC_a36', 'NUM_CONC_a37', 'NUM_CONC_a38', 'NUM_CONC_a39', 'NUM_CONC_a40', 'BIN_CENTERS', 
                'BIN_EDGES','TOT_MASS_CONC', 'TOT_NUM_CONC', 'TOT_WET_NUM_CONC', 'TOT_HYDROPHOBIC_MASS_CONC', 
                'TOT_HYDROPHYLIC_MASS_CONC', 'PM1_MASS_CONC', 'PM25_MASS_CONC', 'PM10_MASS_CONC', 'EXT_AER_550', 
                'EXT_AER_550_INTERNAL', 'EXT_AER_550_EXTERNAL', 'SCAT_AER_550', 'SCAT_AER_550_INTERNAL', 
                'SCAT_AER_550_EXTERNAL', 'NUM_CONC_A1', 'NUM_CONC_A2', 'NUM_CONC_A3', 'MASS_CONC_A1', 'MASS_CONC_A2', 
                'MASS_CONC_A3', 'SCAT_AER_550_PR_A1', 'SCAT_AER_550_INTERNAL_A1', 'SCAT_AER_550_EXTERNAL_A1', 
                'SCAT_AER_550_PR_A2', 'SCAT_AER_550_INTERNAL_A2', 'SCAT_AER_550_EXTERNAL_A2', 'SCAT_AER_550_PR_A3', 
                'SCAT_AER_550_INTERNAL_A3', 'SCAT_AER_550_EXTERNAL_A3', 'EXT_AER_550_PR_A1', 'EXT_AER_550_INTERNAL_A1', 
                'EXT_AER_550_EXTERNAL_A1', 'EXT_AER_550_PR_A2', 'EXT_AER_550_INTERNAL_A2', 'EXT_AER_550_EXTERNAL_A2', 
                'EXT_AER_550_PR_A3', 'EXT_AER_550_INTERNAL_A3', 'EXT_AER_550_EXTERNAL_A3', 'ccn_pr_001_a1', 
                'ccn_pr_001_a2', 'ccn_pr_001_a3', 'ccn_pr_003_a1', 'ccn_pr_003_a2', 'ccn_pr_003_a3', 'ccn_pr_006_a1', 
                'ccn_pr_006_a2', 'ccn_pr_006_a3', 'ccn_pr_010_a1', 'ccn_pr_010_a2', 'ccn_pr_010_a3', 'ccn_internal_001_a1', 
                'ccn_internal_001_a2', 'ccn_internal_001_a3', 'ccn_internal_003_a1', 'ccn_internal_003_a2', 
                'ccn_internal_003_a3', 'ccn_internal_006_a1', 'ccn_internal_006_a2', 'ccn_internal_006_a3', 
                'ccn_internal_010_a1', 'ccn_internal_010_a2', 'ccn_internal_010_a3', 'ccn_external_001_a1', 
                'ccn_external_001_a2', 'ccn_external_001_a3', 'ccn_external_003_a1', 'ccn_external_003_a2', 
                'ccn_external_003_a3', 'ccn_external_006_a1', 'ccn_external_006_a2', 'ccn_external_006_a3', 
                'ccn_external_010_a1', 'ccn_external_010_a2', 'ccn_external_010_a3', 'N_PARTS', 'CELL_VOL', 
                'N_COMPONENTS', 'TOT_NUM_CONC_COAGULATED', 'TOT_BC_NUM_CONC', 'TOT_BC_NUM_CONC_AGED', 
                'TOT_COAGULATION_NUM_CONC', 'D_ALPHA', 'D_GAMMA', 'CHI', 'D_ALPHA_CCN', 'D_GAMMA_CCN', 'CHI_CCN', 
                'D_ALPHA_OPT', 'D_GAMMA_OPT', 'CHI_OPT', 'D_ALPHA_SUBMICRON', 'D_GAMMA_SUBMICRON', 'CHI_SUBMICRON', 
                'D_ALPHA_CCN_SUBMICRON', 'D_GAMMA_CCN_SUBMICRON', 'CHI_CCN_SUBMICRON', 'D_ALPHA_OPT_SUBMICRON', 
                'D_GAMMA_OPT_SUBMICRON', 'CHI_OPT_SUBMICRON', 'D_ALPHA_SPECIES_A1', 'D_GAMMA_SPECIES_A1', 'CHI_SPECIES_A1', 
                'D_ALPHA_CCN_A1', 'D_GAMMA_CCN_A1', 'CHI_CCN_A1', 'D_ALPHA_OPT_A1', 'D_GAMMA_OPT_A1', 'CHI_OPT_A1', 
                'D_ALPHA_SPECIES_A2', 'D_GAMMA_SPECIES_A2', 'CHI_SPECIES_A2', 'D_ALPHA_CCN_A2', 'D_GAMMA_CCN_A2', 
                'CHI_CCN_A2', 'D_ALPHA_OPT_A2', 'D_GAMMA_OPT_A2', 'CHI_OPT_A2', 'D_ALPHA_SPECIES_A3', 'D_GAMMA_SPECIES_A3', 
                'CHI_SPECIES_A3', 'D_ALPHA_CCN_A3', 'D_GAMMA_CCN_A3', 'CHI_CCN_A3', 'D_ALPHA_OPT_A3', 'D_GAMMA_OPT_A3', 
                'CHI_OPT_A3', 'pmc_SO4', 'pmc_NO3', 'pmc_Cl', 'pmc_NH4', 'pmc_MSA', 'pmc_ARO1', 'pmc_ARO2', 'pmc_ALK1', 
                'pmc_OLE1', 'pmc_API1', 'pmc_API2', 'pmc_LIM1', 'pmc_LIM2', 'pmc_CO3', 'pmc_Na', 'pmc_Ca', 'pmc_OIN', 
                'pmc_OC', 'pmc_BC', 'pmc_H2O', 'ccn_001', 'ccn_003', 'ccn_006', 'ccn_010', 'ccn_internal_001', 
                'ccn_internal_003', 'ccn_internal_006', 'ccn_internal_010', 'ccn_external_001', 'ccn_external_003', 
                'ccn_external_006', 'ccn_external_010', 'num_conc_source_000', 'num_conc_source_001', 'num_conc_source_002', 
                'num_conc_source_003']
    gas_vars = ['h2so4', 'hno3', 'hcl', 'nh3', 'no', 'no2', 'no3', 'n2o5', 'hono', 'hno4', 'o3', 'o1d', 'O3P', 'oh', 'ho2', 
                'h2o2', 'co', 'so2', 'ch4', 'c2h6', 'ch3o2', 'ethp', 'hcho', 'ch3oh', 'ANOL', 'ch3ooh', 'ETHOOH', 'ald2', 
                'hcooh', 'RCOOH', 'c2o3', 'pan', 'aro1', 'aro2', 'alk1', 'ole1', 'api1', 'api2', 'lim1', 'lim2', 'par', 
                'AONE', 'mgly', 'eth', 'OLET', 'OLEI', 'tol', 'xyl', 'cres', 'to2', 'cro', 'open', 'onit', 'rooh', 'ro2', 
                'ano2', 'nap', 'xo2', 'xpar', 'isop', 'isoprd', 'isopp', 'isopn', 'isopo2', 'api', 'lim', 'dms', 'msa', 
                'dmso', 'dmso2', 'ch3so2h', 'ch3sch2oo', 'ch3so2', 'ch3so3', 'ch3so2oo', 'ch3so2ch2oo', 'SULFHOX',]

    archive_path = None
    aero_data = {}
    scenario_colors = {}
    aerodist_data = {}
    wrf_data = {}
    scenario_slurm_map = {}
    nsh_dict = {}
    boxplot_data = {}
    gas_fmt_map = {}
    aerosol_fmt_map = {}
    n_times = None
    n_levels = None 
    domain_x_cells = None
    domain_y_cells = None
    historydelta_m = None
    #gridsize = None

    def __init__(self):
        self._getFormatGasSpecies()
        self._getFormatAerosolSpecies()


    def addScenario(self, scenario_name, slurm_id, **kwargs):
        start = kwargs.get('start_time', '2023-03-20_09:00:00')
        self.aero_data[scenario_name] =  nc.Dataset(f'{self.archive_path}/slurm-{slurm_id}/aerosols_d01_{start}')
        self.aerodist_data[scenario_name] =  nc.Dataset(f'{self.archive_path}/slurm-{slurm_id}/aerosol_dist_d01_{start}')
        self.wrf_data[scenario_name] = nc.Dataset(f'{self.archive_path}/slurm-{slurm_id}/wrfout_d01_{start}')
        self.scenario_slurm_map[scenario_name] = slurm_id
        self.nsh_dict[scenario_name] = {}

        if (scenario_name == 'uniform-basecase' and not self.n_times):
            self.addSimAttributes()
    
    def addSimAttributes(self):
        try:
            basecase_aerodata = self.aero_data['uniform-basecase']
            self.n_times = basecase_aerodata.dimensions['Time'].size
            self.n_levels = basecase_aerodata.dimensions['bottom_top'].size
            self.domain_x_cells = basecase_aerodata.dimensions['west_east'].size
            self.domain_y_cells = basecase_aerodata.dimensions['south_north'].size

            times = self.wrf_data['uniform-basecase']['Times'][:].data
            timestamps = [''.join(times[i].astype(str)) for i in range(times.shape[0])]
            timestamps_dt = pd.to_datetime(timestamps, format='%Y-%m-%d_%H:%M:%S')
            self.historydelta_m = (timestamps_dt[1] - timestamps_dt[0]).total_seconds()/60.0

        except KeyError:
            print('uniform-basecase not found in aero_data')
    
    def getScenarioList(self):
        return list(self.aero_data.keys())
    
    def getScenarioSH(self, return_scaling=False):
        path = '/data/keeling/a/sf20/b/wrf-partmc-spatial-het/WRFV3/test/em_les/spatial-het'
        filename = f'sh_patterns_xres{self.domain_x_cells}_yres{self.domain_y_cells}_exact.csv'
        dataset = pd.read_csv(os.path.join(path, filename), index_col='scenario')
        scenario_list = self.getScenarioList()
        if not return_scaling:
            scenario_sh = {}
            for scenario in scenario_list:
                if 'no-nh4' in scenario:
                    proxy_scenario = scenario.replace('-no-nh4', '')
                    scenario_sh[scenario] = dataset.loc[proxy_scenario, 'NSH']
                else:
                    scenario_sh[scenario] = dataset.loc[scenario, 'NSH']
            # sort by ascending SH (or scaling value)
            scenario_sh = {k:v for k, v in sorted(scenario_sh.items(), key=lambda item: item[1])}
            return scenario_sh
        scenario_sh = {}
        scaling_sh = {}
        for scenario in scenario_list:
            if 'no-nh4' in scenario:
                proxy_scenario = scenario.replace('-no-nh4', '')
                scenario_sh[scenario] = dataset.loc[proxy_scenario, 'NSH']
                scaling_sh[scenario] = dataset.loc[proxy_scenario, 'scaling-factor']
            else:
                scenario_sh[scenario] = dataset.loc[scenario, 'NSH']
                scaling_sh[scenario] = dataset.loc[scenario, 'scaling-factor']
            
        # sort by ascending SH (or scaling value)
        scenario_sh = {k:v for k, v in sorted(scenario_sh.items(), key=lambda item: item[1])}
        scaling_sh = {k:v for k, v in sorted(scaling_sh.items(), key=lambda item: item[1])}
        return scenario_sh, scaling_sh
    
    def _formatSpeciesSubscripts(self, var):
        ion_charges = {'NH4': '+', 'NO3': '-', 'SO4': '2-'}
        set_ion_charge = False
        if var in ion_charges:
            set_ion_charge = True
        l = list(var.upper())
        l_orig = l.copy()
        chars_modified = 0
        for i, char in enumerate(l_orig):
            #print(char)
            try:
                int(char)
            except ValueError:
                continue
            #print(i)
            l.insert(i+3*chars_modified, '$')
            #print(l)
            l.insert(i+1+ 3*chars_modified, '_')
            #print(l)
            l.insert(i+3 + 3*chars_modified, '$')
            #print(l)
            chars_modified += 1
        formatted_var = ''.join(l)
        if set_ion_charge:
            formatted_var = formatted_var[:-1] + '^{' + ion_charges[var] + '}$'
        l = []
        return formatted_var

    def _getFormatGasSpecies(self):
        for var in self.gas_vars:
            self.gas_fmt_map[var] = self._formatSpeciesSubscripts(var)

    def _getFormatAerosolSpecies(self):
        ccn_vars = {'D_ALPHA': '$D_{\\alpha}$', 
                    'D_GAMMA': '$D_{\gamma}$',  
                    'CHI': '$\chi$', 
                    'D_ALPHA_CCN': '$D_{\\alpha\mathrm{, CCN}}$',  
                    'D_GAMMA_CCN': '$D_{\gamma\mathrm{, CCN}}$',  
                    'CHI_CCN': '$\chi_{\mathrm{CCN}}$'
                    }
        for var in self.aero_vars:
            if var.startswith('pmc_'):
                var_temp = var.replace('pmc_', '')
                fmt_var = self._formatSpeciesSubscripts(var_temp)
                self.aerosol_fmt_map[var] = ''.join(['Aerosol ', fmt_var])
            elif var.startswith('ccn'):
                var_temp = var.upper().split('_')
                if len(var_temp) == 2:
                    supersat = int(var_temp[-1]) / 10
                    var_fmt = f'{var_temp[0]} ($S={supersat}\%$)'
                    self.aerosol_fmt_map[var] = var_fmt
                elif len(var_temp) == 3:
                    supersat = int(var_temp[-1]) / 10
                    var_fmt = f'{var_temp[0]} {var_temp[1].title()} ($S={supersat}\%$)'
                    self.aerosol_fmt_map[var] = var_fmt
                else:
                   self.aerosol_fmt_map[var] = var
            elif var in ccn_vars:
                self.aerosol_fmt_map[var] = ccn_vars[var]
            else:
                self.aerosol_fmt_map[var] = var
    
    def getScenarioGeneralLabels(self):
        scenario_sh = self.getScenarioSH()
        scenario_labels = {}
        i = 0
        for scenario in scenario_sh.keys():
            if scenario == 'uniform-basecase':
                scenario_labels[scenario] = 'Uniform base case'
                i += 1
            elif 'no-nh4' in scenario:
                base_scenario = scenario.replace('-no-nh4', '')
                try:
                    base_scenario_general_label = scenario_labels[base_scenario]
                except KeyError:
                    print(f"No scenario {base_scenario} found - did you load the nh4 scenario and not the run with nh4?")
                scenario_labels[scenario] = f'{base_scenario_general_label}, no NH$_4$'
            else:
                scenario_labels[scenario] = f'Scenario {i}'
                i += 1
        return scenario_labels
    
    def getScenarioColors(self):
        
        scenario_colors_standard = {"uniform-basecase": "",
                           "fx2fy2": "", 
                           "fx1fy0": "", 
                           "road-10x": "", 
                           "point-source-10x10": "", 
                           "point-source-1x1": "",
                           }
        scenario_list = list(scenario_colors_standard.keys())
        viridis = plt.get_cmap('viridis')
        color_vals = np.linspace(0.1, .9, len(scenario_colors_standard.keys())) # just a fix
        
        # kind of hacky but working solution
        scenario_colors = {"uniform-basecase": "",
                           "uniform-basecase-no-nh4": "",
                           "fx2fy2": "", 
                           "fx1fy0": "", 
                           "road-10x": "", 
                           "point-source-10x10": "", 
                           "point-source-1x1": "",
                           "point-source-1x1-no-nh4": ""
                           }
        
        for i, scenario in zip(color_vals, scenario_colors_standard.keys()):
            if scenario == 'uniform-basecase':
                scenario_colors[scenario] = 'k'
                scenario_colors['uniform-basecase-no-nh4'] = 'k'
            else:
                rgba = viridis(i)
                hex_color = mplcolors.to_hex(rgba)
                scenario_colors[scenario] = hex_color
                if scenario == 'point-source-1x1':
                    rgba = viridis(i)
                    hex_color = mplcolors.to_hex(rgba)
                    scenario_colors['point-source-1x1-no-nh4'] = hex_color

        return scenario_colors


global Archive
Archive = DataStruct()


class GriddedOutputDataStruct(DataStruct):
    bin_edges = None
    bin_logwidth = None
    bin_geocenter = None
    n_bins = None
    bin_min = None
    bin_max = None
    gridded_data = {}

    #sim_dict = None  # replace with wrf_data?
    
    #t_idx = None # NOTE: what is this for?

    # Try using the inherited dictionaries from the DataStruct class
    #gas_species = None
    #aero_species = None

    # use domain_x_cells and domain_y_cells from inherited DataStruct class
    #nx = None
    #ny = None

    def __init__(self, **kwargs):
        # Assume default WRF-PartMC binning scheme (100 bins ranging from 1e-9 to 1e-3 m)
        self.bin_min= -9
        self.bin_max = -3
        self.n_bins = 101
        self.bin_edges = np.logspace(self.bin_min, self.bin_max, self.n_bins)
        self.bin_logwidth = (np.log10(self.bin_edges[1:]) - np.log10(self.bin_edges[0:-1]))[0]
        self.bin_geocenter = np.sqrt(self.bin_edges[:-1]*self.bin_edges[1:])
        
        # if necessary
        self.aero_species = [var.replace('pmc_', '').upper() for var in self.aero_vars if var.startswith('pmc')]
        self.gas_species = [var.upper() for var in self.gas_vars]
    
    def _getVertGridCellPartIndices(self, data, k):
        # return the start and end indices of the particles in the requested kth vertical 
        # grid cell
        # TODO: with this method I'm getting one less particle in the topmost grid cell
        # when I check against data['n_parts'][-1].item()
        #k -= 1 # convert from one indexing to zero indexing 

        start_idx = data['part_start_index'][k]
        n_parts_in_cell = data['n_parts'][k].item()
        end_idx = start_idx + n_parts_in_cell
        return start_idx, end_idx
    
    def _getParticleDiameters(self, data):
        particle_volume = (data['aero_particle_mass'][:].T/data['aero_density'][:].reshape(1, 20)).sum(axis=1)
        particle_diameter = np.cbrt(particle_volume*(6/np.pi))
        return particle_diameter


    def _computeParticleKappa(self, data):
        aero_density = data['aero_density']
        aero_kappa = data['aero_kappa']

        # volume fraction of each solid species (every species except last component which is water)
        aero_solid_volume = (data['aero_particle_mass'][:] /aero_density[:].data[:, np.newaxis]).data[:-1, :]

        # sum up the volumes for all solid species 
        aero_solid_total_volume = aero_solid_volume.sum(axis=0)

        # volume fraction of each solid species
        aero_solid_volume_frac = aero_solid_volume/aero_solid_total_volume

        kappa = (aero_solid_volume_frac * aero_kappa[:].data[:-1, np.newaxis]).sum(axis=0)

        return kappa

    def processCrossSection(self, data_path, xstart, xend, ystart, yend, z_idx, t_idx):

        xwidth = xend-xstart + 1
        ywidth = yend-ystart + 1

        crosssec_aero_diams = np.array([])
        crosssec_aero_masses = np.zeros((20, 1))
        crosssec_aero_numconc = np.array([])
        crosssec_aero_kappa = np.array([])

        crosssec_aero_component_start_ind = np.zeros((xwidth, ywidth)).astype(int)
        cell_start_idx = 0
        for x_idx in np.arange(xwidth):
            for y_idx in np.arange(ywidth):
                x_cell = xstart + x_idx
                y_cell = ystart + y_idx
                filename = f'gridded-output_{str(x_cell).zfill(3)}_{str(y_cell).zfill(3)}_{str(t_idx).zfill(8)}.nc'

                data = nc.Dataset(os.path.join(data_path, filename))

                # Retreive aerosol particle array indices for vertical level
                start_idx, end_idx = self._getVertGridCellPartIndices(data, z_idx)

                # Species attributes
                # Aerosol species density
                species_density = data['aero_density'][:]
                # Aerosol species kappa
                species_kappa = data['aero_kappa'][:]
                # Aerosol species molecular weight
                species_molec_weight = data['aero_molec_weight'][:]
                # Number of ions after dissocation of aerosol species
                species_num_ions = data['aero_num_ions'][:]

                # Aerosol Diameters
                particle_diameters = self._getParticleDiameters(data)
                level_part_diams = particle_diameters[start_idx:end_idx]
                #print(cell_part_diams[:1])
                # Store the aerosol particle diameters in a 1D array for the specified vertical level
                crosssec_aero_diams = np.append(crosssec_aero_diams, level_part_diams)

                # Aerosol Masses
                level_part_masses = data['aero_particle_mass'][:, start_idx:end_idx] # constituent masses of each aerosol particle
                #print(cell_part_masses[:, :1][0])
                # Store the aerosol particle masses in a 20 X N for the specified vertical level (N total number of aerosol particles in the level)
                # 20 constituent species masses in each aerosol particle
                crosssec_aero_masses = np.append(crosssec_aero_masses, level_part_masses, axis=1)

                # Aerosol number concentration
                level_aero_numconc = data['aero_num_conc'][start_idx:end_idx]
                crosssec_aero_numconc = np.append(crosssec_aero_numconc, level_aero_numconc)
                
                # Aerosol kappa (hygroscopicity)
                aero_kappa = self._computeParticleKappa(data)
                cell_aero_kappa = aero_kappa[start_idx:end_idx]
                crosssec_aero_kappa = np.append(crosssec_aero_kappa, cell_aero_kappa)
                
                # Store location to the starting index for the aerosol particle in each grid cell (dimensions nx x ny for accessing aerosol particles in each cell)
                crosssec_aero_component_start_ind[x_idx, y_idx] = cell_start_idx
                cell_start_idx  = cell_start_idx + (end_idx-start_idx)


        crosssec_aero_masses = crosssec_aero_masses[:, 1:] # need to remove the first element since array of zeros used as placeholder 

        return (crosssec_aero_diams, crosssec_aero_numconc, crosssec_aero_masses, crosssec_aero_kappa, species_density,
                species_kappa, species_molec_weight, species_num_ions)
    
    def loadData(self, scenario, xstart, xend, ystart, yend, z_idx, t_idx, verbose=True):
        slurmid = self.scenario_slurm_map[scenario]
        output_path = self.archive_path
        output_subdir = os.path.join(output_path, f'slurm-{slurmid}')

        if os.path.isfile(f'./crosssec_data/crosssec_{scenario}_t{t_idx}_z{z_idx}.nc'):
            if verbose:
                print('Loading file')
            crosssec_data = nc.Dataset(f'./crosssec_data/crosssec_{scenario}_t{t_idx}_z{z_idx}.nc', 'r', format='NETCDF4')
            # load the environmental variables
            crosssec_aero_diams = crosssec_data['aero_diams'][:]
            crosssec_aero_numconc = crosssec_data['aero_numconc'][:]
            crosssec_aero_masses = crosssec_data['aero_masses'][:]
            crosssec_aero_kappa = crosssec_data['aero_kappa'][:]
            crosssec_species_density = crosssec_data['species_density'][:]
            crosssec_species_kappa = crosssec_data['species_kappa'][:]
            crosssec_species_molec_weight = crosssec_data['species_molec_weight'][:]
            crosssec_species_num_ions = crosssec_data['species_num_ions'][:] 

            crosssec_data.close()
        else:
            if verbose:
                print('File does not exist, processing data')
            data_tuple = self.processCrossSection(output_subdir, xstart, xend, ystart, yend, z_idx, t_idx)
            (crosssec_aero_diams, crosssec_aero_numconc, crosssec_aero_masses, crosssec_aero_kappa, crosssec_species_density,
            crosssec_species_kappa, crosssec_species_molec_weight, crosssec_species_num_ions) = data_tuple

            history_dt = self.historydelta_m/60 # hours
            time =  (t_idx-1)*history_dt
            processed_data = nc.Dataset(f'./crosssec_data/crosssec_{scenario}_t{t_idx}_z{z_idx}.nc', 'w', format='NETCDF4')
            processed_data.description = f'Processed cross-section simulation data at time = {time} hr, zlevel = {z_idx}'
            # dimensions
            #time_dim_size = 7
            #processed_data.createDimension('time', time_dim_size)
            n_particles = crosssec_aero_diams.size
            processed_data.createDimension('n_particles', n_particles) #NOTE the issue with this approach is that the number of particle changes with time.
            n_species = len(self.aero_species)
            processed_data.createDimension('n_species', n_species)
            aero_diams = processed_data.createVariable('aero_diams', 'f8', ('n_particles'))
            aero_numconc = processed_data.createVariable('aero_numconc', 'f8', ('n_particles'))
            aero_masses = processed_data.createVariable('aero_masses', 'f8', ('n_species', 'n_particles'))
            aero_kappa = processed_data.createVariable('aero_kappa', 'f8', ('n_particles'))
            species_density = processed_data.createVariable('species_density', 'f8', ('n_species'))
            species_kappa = processed_data.createVariable('species_kappa', 'f8', ('n_species'))
            species_molec_weight = processed_data.createVariable('species_molec_weight', 'f8', ('n_species'))
            species_num_ions = processed_data.createVariable('species_num_ions', 'f8', ('n_species'))
            
            # variables
            aero_diams[:] = crosssec_aero_diams
            aero_numconc[:] = crosssec_aero_numconc
            aero_masses[:] = crosssec_aero_masses
            aero_kappa[:] = crosssec_aero_kappa
            species_density[:] = crosssec_species_density
            species_kappa[:] = crosssec_species_kappa
            species_molec_weight[:] = crosssec_species_molec_weight
            species_num_ions[:] = crosssec_species_num_ions


            processed_data.close()
        
        xwidth = xend-xstart + 1
        ywidth = yend-ystart + 1
        n_total_cells = xwidth*ywidth

        self.gridded_data[scenario] = {}
        self.gridded_data[scenario]['xstart'] = xstart
        self.gridded_data[scenario]['xend'] = xend
        self.gridded_data[scenario]['ystart'] = ystart
        self.gridded_data[scenario]['yend'] = yend
        self.gridded_data[scenario]['n_total_cells'] = n_total_cells
        self.gridded_data[scenario]['aero_diams'] = crosssec_aero_diams
        self.gridded_data[scenario]['aero_numconc'] = crosssec_aero_numconc
        self.gridded_data[scenario]['aero_masses'] = crosssec_aero_masses
        self.gridded_data[scenario]['aero_kappa'] = crosssec_aero_kappa
        self.gridded_data[scenario]['species_density'] = crosssec_species_density
        self.gridded_data[scenario]['species_kappa'] = crosssec_species_kappa
        self.gridded_data[scenario]['species_molec_weight'] = crosssec_species_molec_weight
        self.gridded_data[scenario]['species_num_ions'] = crosssec_species_num_ions
        #return crosssec_aero_diams, crosssec_aero_numconc, crosssec_aero_masses

    
global GriddedOutput
GriddedOutput = GriddedOutputDataStruct()