import os
import netCDF4 as nc

output_path = '/data/nriemer/d/sf20/les_output/wrf-partmc'

# Uniform basecase 
print('..loading basecase')
basecase_subdir = os.path.join(output_path, 'slurm-1950592')
basecase_aerodata = nc.Dataset(os.path.join(basecase_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
#basecase_aerodistdata = nc.Dataset(os.path.join(basecase_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

#point-source-10x10 (tchem = 1 s)
print('..loading point source')
pointsource10x10_subdir = os.path.join(output_path, 'slurm-1951163')
pointsource10x10_aerodata = nc.Dataset(os.path.join(pointsource10x10_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
#pointsource10x10_aerodistdata = nc.Dataset(os.path.join(pointsource10x10_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

#point-source-10x10, tchem = 60 seconds 
#print('..loading point source')
pointsource10x10tchem60_subdir = os.path.join(output_path, 'slurm-1953691')
pointsource10x10tchem60_aerodata = nc.Dataset(os.path.join(pointsource10x10tchem60_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
#pointsource10x10tchem60_aerodistdata = nc.Dataset(os.path.join(pointsource10x10tchem60_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))


global aerodata_dict
global nsh_dict

aerodata_dict = {'basecase': {'aerodata': basecase_aerodata,
                              #'distdata': basecase_aerodistdata,
                              },
                 'point-source-10x10': {'aerodata': pointsource10x10_aerodata,
                                 #'distdata': pointsource10x10_aerodistdata,
                                 },
                  'point-source-10x10-tchem60': {'aerodata': pointsource10x10tchem60_aerodata,
                                 #'distdata': pointsource10x10tchem60_aerodistdata,
                                 },
                }

nsh_dict = {'basecase': {},
            'point-source-10x10': {},
            'point-source-10x10-tchem60': {},
            }

n_times = basecase_aerodata.dimensions['Time'].size
n_levels = basecase_aerodata.dimensions['bottom_top'].size
domain_x_cells = basecase_aerodata.dimensions['west_east'].size
domain_y_cells = basecase_aerodata.dimensions['south_north'].size

global boxplot_data
boxplot_data = {}

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
   
print('Finished loading data structures')