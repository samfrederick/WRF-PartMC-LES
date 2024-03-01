import netCDF4 as nc
import numpy as np
import os
from sys import getsizeof

global bin_edges, bin_logwidth, n_bins, gas_species, aero_species, sim_dict

nx, ny = 40, 40
t_idx = 1

sim_dict = {'basecase': 1909559, # aerosol emissions begin after 1 hr spinup
            #'basecase': 1926393, # uniform gas IC matches aerosol IC, aerosol emissions begin after 1 hr spinup
            #'point-source-10x10': 1908083 #  aerosol emissions begin at sim start
            'point-source-10x10': 1919124  # aerosol emissions begin after 1 hr spinup #1909980 INCORRECT GAS EMISSION RATES
            #'point-source-10x10':1927316 # uniform gas IC matches aerosol IC, aerosol emissions begin after 1 hr spinup
            }

BIN_SCHEME = 'HI_RES'
#BIN_SCHEME = 'WRF_PARTMC'

if BIN_SCHEME == 'HI_RES':
    bin_min= -8
    bin_max = -5
    n_bins = 101

# Matching default WRF-PartMC binning scheme
if BIN_SCHEME == 'WRF_PARTMC':
    bin_min= -9
    bin_max = -3
    n_bins = 101

bin_edges = np.logspace(bin_min, bin_max, n_bins)
bin_logwidth = (np.log10(bin_edges[1:]) - np.log10(bin_edges[0:-1]))[0]
bin_geocenter = np.sqrt(bin_edges[:-1]*bin_edges[1:])

gas_species = ['H2SO4',
                'HNO3',
                'HCl',
                'NH3',
                'NO',
                'NO2',
                'NO3',
                'N2O5',
                'HONO',
                'HNO4',
                'O3',
                'O1D',
                'O3P',
                'OH',
                'HO2',
                'H2O2',
                'CO',
                'SO2',
                'CH4',
                'C2H6',
                'CH3O2',
                'ETHP',
                'HCHO',
                'CH3OH',
                'ANOL',
                'CH3OOH',
                'ETHOOH',
                'ALD2',
                'HCOOH',
                'RCOOH',
                'C2O3',
                'PAN',
                'ARO1',
                'ARO2',
                'ALK1',
                'OLE1',
                'API1',
                'API2',
                'LIM1',
                'LIM2',
                'PAR',
                'AONE',
                'MGLY',
                'ETH',
                'OLET',
                'OLEI',
                'TOL',
                'XYL',
                'CRES',
                'TO2',
                'CRO',
                'OPEN',
                'ONIT',
                'ROOH',
                'RO2',
                'ANO2',
                'NAP',
                'XO2',
                'XPAR',
                'ISOP',
                'ISOPRD',
                'ISOPP',
                'ISOPN',
                'ISOPO2',
                'API',
                'LIM',
                'DMS',
                'MSA',
                'DMSO',
                'DMSO2',
                'CH3SO2H',
                'CH3SCH2OO',
                'CH3SO2',
                'CH3SO3',
                'CH3SO2OO',
                'CH3SO2CH2OO',
                'SULFHOX']
aero_species = ['SO4',
                'NO3',
                'Cl',
                'NH4',
                'MSA',
                'ARO1',
                'ARO2',
                'ALK1',
                'OLE1',
                'API1',
                'API2',
                'LIM1',
                'LIM2',
                'CO3',
                'Na',
                'Ca',
                'OIN',
                'OC',
                'BC',
                'H2O']


"""
# load coarse-resolved datasets
output_path = '/data/nriemer/d/sf20/les_output/wrf-partmc'
basecase_subdir = os.path.join(output_path, 'slurm-1851783')#'slurm-1826399') 
basecase_aerodata = nc.Dataset(os.path.join(basecase_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
basecase_aerodistdata = nc.Dataset(os.path.join(basecase_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

pointsource10x10_subdir = os.path.join(output_path, 'slurm-1853275') # aerosol ICs everywhere, emissions follow SH patttern
pointsource10x10_aerodata = nc.Dataset(os.path.join(pointsource10x10_subdir, 'aerosols_d01_2023-06-21_09:00:00'))
pointsource10x10_aerodistdata = nc.Dataset(os.path.join(pointsource10x10_subdir, 'aerosol_dist_d01_2023-06-21_09:00:00'))

global aerodata_dict
global nsh_dict

aerodata_dict = {
                 
                 'basecase': {'aerodata': basecase_aerodata,
                              'distdata': basecase_aerodistdata,
                              },

                 'point-source-10x10': {'aerodata': pointsource10x10_aerodata,
                                 'distdata': pointsource10x10_aerodistdata,
                                 },
                }

nsh_dict = {'basecase': {},
            
            'point-source-10x10': {},

            }

n_times = basecase_aerodata.dimensions['Time'].size
n_levels = basecase_aerodata.dimensions['bottom_top'].size
domain_x_cells = basecase_aerodata.dimensions['west_east'].size
domain_y_cells = basecase_aerodata.dimensions['south_north'].size
"""
