## Repository for WRF-PartMC-LES model initialization code and analysis

This repository contains scipts and documentation for the use of WRF-PartMC configured
for large-eddy simulations (LES). The meteorological initial conditions and relevant namelist
parameters for LES are adopted from the Weather Research and Forecasting (WRF) model's LES 
test case. 

### Recommended installation procedure

WRF and the Particle Monte Carlo (PartMC) aerosol model are coupled specifically for WRF v3.9.1. Please
see the [WRF-PartMC repository](https://github.com/open-atmos/wrf-partmc/tree/v1.0) for details on its installation 
and associated dependencies first (including PartMC and the multiphase chemistry module MOSAIC). WRF-PartMC 
should be compiled for the LES test case prior to running scripts in this repository (namely, the `wrf-partmc_run.sh`
bash script for submitting WRF-PartMC-LES jobs to SLURM which requires the `wrf.exe` and `ideal.exe` executables to be
present in the main directory of this repository following compilation of WRF-PartMC). 

After cloning the WRF-PartMC repository and related dependencies and compiling the code for the LES test case, navigate 
to the LES test case subdirectory (`~/wrf-partmc/WRFV3/test/em_les`). Ensure that the main WRF-PartMC executables `wrf.exe` and `ideal.exe`
are present. This repository should be cloned at this subdirectory path. 