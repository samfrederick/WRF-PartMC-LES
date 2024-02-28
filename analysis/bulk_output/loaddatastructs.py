import os
import netCDF4 as nc

output_path = '/data/nriemer/d/sf20/les_output/wrf-partmc'
# 1826399 basecase without chemistry
# 1851783 basecase with mosaic chem, tstart 2023-06-21_09:00:00
# 1909559, basecase aerosol emissions delayed by 1 hour, tstart 2023-06-21_09:00:00
# 1934602 basecase, corrected IC aitken mode geom mean diam, corrected IC gas concentrations, tstart 2023-03-20_09:00:00
basecase_subdir = os.path.join(output_path, 'slurm-1934602')
basecase_aerodata = nc.Dataset(os.path.join(basecase_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
basecase_aerodistdata = nc.Dataset(os.path.join(basecase_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

"""
# Basecase aerosol emissions delayed by 1 hour, change date to vernal equinox
basecasedelayemiss_subdir = os.path.join(output_path, 'slurm-1909559')
basecasedelayemiss_aerodata = nc.Dataset(os.path.join(basecasedelayemiss_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
basecasedelayemiss_aerodistdata = nc.Dataset(os.path.join(basecasedelayemiss_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

# Basecase NO CHEM
basecasenochem_subdir = os.path.join(output_path, 'slurm-1826399')
basecasenochem_aerodata = nc.Dataset(os.path.join(basecasenochem_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
basecasenochem_aerodistdata = nc.Dataset(os.path.join(basecasenochem_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

# Basecase NO COAGULATION
basecasenocoag_subdir = os.path.join(output_path, 'slurm-1901327')
basecasenocoag_aerodata = nc.Dataset(os.path.join(basecasenocoag_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
basecasenocoag_aerodistdata = nc.Dataset(os.path.join(basecasenocoag_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))
"""

# aerosol ICS follow emission ICS: 1827809
# aerosol ICs everywhere, no chem: 1835024
# scenario with mosiac chem 1852006
fx1fy0_subdir = os.path.join(output_path, 'slurm-1852006')
fx1fy0_aerodata = nc.Dataset(os.path.join(fx1fy0_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
fx1fy0_aerodistdata = nc.Dataset(os.path.join(fx1fy0_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

#pointsource_subdir = os.path.join(output_path, 'slurm-1848762') # aerosol ICs everywhere, emissions follow SH patttern
#pointsource_aerodata = nc.Dataset(os.path.join(pointsource_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
#pointsource_aerodistdata = nc.Dataset(os.path.join(pointsource_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

#point-source-10x10 - 1853275 NOTE SHOULD BE DEPRECATED DUE TO INCORRECT RATE OF GAS EMISSIONS (NOT SCALED), tstart 2023-06-21_09:00:00
# 1935061, corrected gas emission scaling, corrected IC aitken mode geom mean diam, corrected IC gas concentrations, tstart 2023-03-20_09:00:00
pointsource10x10_subdir = os.path.join(output_path, 'slurm-1935061') # aerosol ICs everywhere, emissions follow SH patttern
pointsource10x10_aerodata = nc.Dataset(os.path.join(pointsource10x10_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
pointsource10x10_aerodistdata = nc.Dataset(os.path.join(pointsource10x10_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

"""
#point-source-10x10, aerosol emissions delayed by 1 hour, change date to vernal equinox
pointsource10x10delayemiss_subdir = os.path.join(output_path, 'slurm-1909980')
pointsource10x10delayemiss_aerodata = nc.Dataset(os.path.join(pointsource10x10delayemiss_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
pointsource10x10delayemiss_aerodistdata = nc.Dataset(os.path.join(pointsource10x10delayemiss_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

#point-source-10x10, aerosol emissions delayed by 1 hour, CORRECTED GAS EMISSIONS
pointsource10x10correctemiss_subdir = os.path.join(output_path, 'slurm-1919124')
pointsource10x10correctemiss_aerodata = nc.Dataset(os.path.join(pointsource10x10correctemiss_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
pointsource10x10correctemiss_aerodistdata = nc.Dataset(os.path.join(pointsource10x10correctemiss_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

#point-source-10x10, NO COAGULATION
pointsource10x10nocoag_subdir = os.path.join(output_path, 'slurm-1896888') # aerosol ICs everywhere, emissions follow SH patttern
pointsource10x10nocoag_aerodata = nc.Dataset(os.path.join(pointsource10x10nocoag_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
pointsource10x10nocoag_aerodistdata = nc.Dataset(os.path.join(pointsource10x10nocoag_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

#point-source-10x10, NO DEPOSITION
pointsource10x10nodep_subdir = os.path.join(output_path, 'slurm-1903184') # aerosol ICs everywhere, emissions follow SH patttern
pointsource10x10nodep_aerodata = nc.Dataset(os.path.join(pointsource10x10nodep_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
pointsource10x10nodep_aerodistdata = nc.Dataset(os.path.join(pointsource10x10nodep_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))
"""

# 1935620, corrected gas emission scaling, corrected IC aitken mode geom mean diam, corrected IC gas concentrations, tstart 2023-03-20_09:00:00
pointsource1x1_subdir = os.path.join(output_path, 'slurm-1935620') # aerosol ICs everywhere, emissions follow SH patttern
pointsource1x1_aerodata = nc.Dataset(os.path.join(pointsource1x1_subdir, 'aerosols_d01_2023-03-20_09:00:00'))
pointsource1x1_aerodistdata = nc.Dataset(os.path.join(pointsource1x1_subdir, 'aerosol_dist_d01_2023-03-20_09:00:00'))

# aerosol ICs everywhere, no chem: 1828109
# scenario with mosiac chem 1852335
fx1fy1_subdir = os.path.join(output_path, 'slurm-1852335')
fx1fy1_aerodata = nc.Dataset(os.path.join(fx1fy1_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
fx1fy1_aerodistdata = nc.Dataset(os.path.join(fx1fy1_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

fx2fy2_subdir = os.path.join(output_path, 'slurm-1854994')
fx2fy2_aerodata = nc.Dataset(os.path.join(fx2fy2_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
fx2fy2_aerodistdata = nc.Dataset(os.path.join(fx2fy2_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

roaddouble_subdir = os.path.join(output_path, 'slurm-1852604')
roaddouble_aerodata = nc.Dataset(os.path.join(roaddouble_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
roaddouble_aerodistdata = nc.Dataset(os.path.join(roaddouble_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

road8x_subdir = os.path.join(output_path, 'slurm-1854165')
road8x_aerodata = nc.Dataset(os.path.join(road8x_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
road8x_aerodistdata = nc.Dataset(os.path.join(road8x_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

road16x_subdir = os.path.join(output_path, 'slurm-1854537')
road16x_aerodata = nc.Dataset(os.path.join(road16x_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
road16x_aerodistdata = nc.Dataset(os.path.join(road16x_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

global aerodata_dict
global nsh_dict

aerodata_dict = {'basecase': {'aerodata': basecase_aerodata,
                              'distdata': basecase_aerodistdata,
                              },
                 #'basecase-delayemiss': {'aerodata': basecasedelayemiss_aerodata,
                 #             'distdata': basecasedelayemiss_aerodistdata,
                 #             },
                 #'basecase-nochem': {'aerodata': basecasenochem_aerodata,
                 #             'distdata': basecasenochem_aerodistdata,
                 #             },
                 #'basecase-nocoag': {'aerodata': basecasenocoag_aerodata,
                 #             'distdata': basecasenocoag_aerodistdata,
                 #             },
                 'fx1fy0': {'aerodata': fx1fy0_aerodata,
                            'distdata': fx1fy0_aerodistdata,
                            },
                 'fx1fy1': {'aerodata': fx1fy1_aerodata,
                            'distdata': fx1fy1_aerodistdata,
                            },
                 'fx2fy2': {'aerodata': fx2fy2_aerodata,
                            'distdata': fx2fy2_aerodistdata,
                            },
                 'road-double': {'aerodata': roaddouble_aerodata,
                            'distdata': roaddouble_aerodistdata,
                            },
                 'road-8x': {'aerodata': road8x_aerodata,
                            'distdata': road8x_aerodistdata,
                            },
                 'road-16x': {'aerodata': road16x_aerodata,
                            'distdata': road16x_aerodistdata,
                            },
                 'point-source-10x10': {'aerodata': pointsource10x10_aerodata,
                                 'distdata': pointsource10x10_aerodistdata,
                                 },
                #'point-source-10x10-delayemiss': {'aerodata': pointsource10x10delayemiss_aerodata,
                 #               'distdata': pointsource10x10delayemiss_aerodistdata,
                 #               },
                 #'point-source-10x10-correctemiss': {'aerodata': pointsource10x10correctemiss_aerodata,
                 #                'distdata': pointsource10x10correctemiss_aerodistdata,
                 #                },
                 #'point-source-10x10-nocoag': {'aerodata': pointsource10x10nocoag_aerodata,
                 #                'distdata': pointsource10x10nocoag_aerodistdata,
                 #                },
                 #'point-source-10x10-nodep': {'aerodata': pointsource10x10nodep_aerodata,
                 #                'distdata': pointsource10x10nodep_aerodistdata,
                 #                },
                 'point-source-1x1': {'aerodata': pointsource1x1_aerodata,
                                 'distdata': pointsource1x1_aerodistdata,
                                 },
                
                 #'fx1fy1': fx1fy1_aerodata
                }

nsh_dict = {'basecase': {},
            #'basecase-delayemiss': {},
            #'basecase-nochem': {},
            #'basecase-nocoag': {},
            'fx1fy0': {},
            'fx1fy1': {},
            'fx2fy2': {},
            'road-double': {},
            'road-8x': {},
            'road-16x': {},
            'point-source-10x10': {},
            #'point-source-10x10-delayemiss': {},
            #'point-source-10x10-correctemiss': {},
            #'point-source-10x10-nocoag': {},
            #'point-source-10x10-nodep': {},
            'point-source-1x1': {},
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
   
  