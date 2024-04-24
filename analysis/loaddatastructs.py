import os
import netCDF4 as nc
import pandas as pd

class DataStruct:
    wrf_vars = ['T', 'P', 'ALT', 'PB', 'DNW', 'DN', 'Z', 'Z_AT_W', 'MAPFAC_M', 'MAPFAC_U', 'MAPFAC_V', 'MAPFAC_MX', 
            'MAPFAC_MY', 'MAPFAC_UX', 'MAPFAC_UY', 'MAPFAC_VX', 'MF_VX_INV', 'MAPFAC_VY','DENSITY_DRY_AIR', 
            'TEMPERATURE', 'REL_HUMID',]
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
    aerodist_data = {}
    wrf_data = {}
    nsh_dict = {}
    boxplot_data = {}
    gas_fmt_map = {}
    aerosol_fmt_map = {}
    n_times = None
    n_levels = None 
    domain_x_cells = None
    domain_y_cells = None
    historydelta_m = None
    gridsize = None

    def __init__(self):
        self._getFormatGasSpecies()
        self._getFormatAerosolSpecies()


    def addScenario(self, scenario_name, slurm_id, **kwargs):
        start = kwargs.get('start_time', '2023-03-20_09:00:00')
        self.aero_data[scenario_name] =  nc.Dataset(f'{self.archive_path}/slurm-{slurm_id}/aerosols_d01_{start}')
        self.aerodist_data[scenario_name] =  nc.Dataset(f'{self.archive_path}/slurm-{slurm_id}/aerosol_dist_d01_{start}')
        self.wrf_data[scenario_name] = nc.Dataset(f'{self.archive_path}/slurm-{slurm_id}/wrfout_d01_{start}')
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
        filename = f'sh_patterns_xres{self.gridsize}_yres{self.gridsize}_exact.csv'
        dataset = pd.read_csv(os.path.join(path, filename), index_col='scenario')
        scenario_list = self.getScenarioList()
        if not return_scaling:
            scenario_sh = {}
            for scenario in scenario_list:
                scenario_sh[scenario] = dataset.loc[scenario, 'NSH']
            # sort by ascending SH (or scaling value)
            scenario_sh = {k:v for k, v in sorted(scenario_sh.items(), key=lambda item: item[1])}
            return scenario_sh
        scenario_sh = {}
        scaling_sh = {}
        for scenario in scenario_list:
            scenario_sh[scenario] = dataset.loc[scenario, 'NSH']
            scaling_sh[scenario] = dataset.loc[scenario, 'scaling-factor']

        # sort by ascending SH (or scaling value)
        scenario_sh = {k:v for k, v in sorted(scenario_sh.items(), key=lambda item: item[1])}
        scaling_sh = {k:v for k, v in sorted(scaling_sh.items(), key=lambda item: item[1])}
        return scenario_sh, scaling_sh
    
    def _formatSpeciesSubscripts(self, var):
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
        l = []
        return formatted_var

    def _getFormatGasSpecies(self):
        for var in self.gas_vars:
            self.gas_fmt_map[var] = self._formatSpeciesSubscripts(var)

    def _getFormatAerosolSpecies(self):
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
            else:
                self.aerosol_fmt_map[var] = var
    
    def getScenarioGeneralLabels(self):
        scenario_sh = self.getScenarioSH()
        scenario_labels = {}
        for i, scenario in enumerate(scenario_sh.keys()):
            scenario_labels[scenario] = f'Scenario {i}'
        return scenario_labels

global Archive
Archive = DataStruct()


