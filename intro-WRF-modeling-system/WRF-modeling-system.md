# What is WRF?
The Weather Research and Forecasting Model (WRF) is a numerical weather prediction system designed to serve both atmospheric research and operational forecasting/hindcasting applications. WRF-Chem is the Weather Research and Forecasting (WRF) model coupled with Chemistry. The model simulates the emission, transport, mixing, and chemical transformation of trace gases and aerosols simultaneously with the meteorology.  

# WRF code distribution
WRF is distrubted at: https://github.com/wrf-model/WRF
The current version of WRF and WRF-Chem are distributed together, if you are running WRF vs. WRF-Chem depends on how you compile the model and your WRF namelist (main input file) and other input files provided.
The current of the WRF modeling system - as of 11 March 2024 is 4.5.2.  

All information on the WRF modleing system including user mailing lists, documentation, and support from NCAR the main model developer for WRF dynamics is found at: https://www2.mmm.ucar.edu/wrf/users/download/get_source.html  
WRF-Chem is managed by NOAA and the main page for information is: https://ruc.noaa.gov/wrf/wrf-chem/

WRF and WRF-Chem are modular modeling systems, which mutliple options for treating processes in the atmosphere, coupling with the land/ocean/ice surfaces, and atmospheric chemistry (trace gas and aerosols).  We will present one model setup that we frequenly use at IGE/LATMOS.  However, other model setups are possible using other atmopserhic, surface, and atmospheric chemsitry setups.

The specific setup for the chemistry part of WRF-Chem presented here is MOZART gas-phase chemistry with MOSAIC aerosols (https://www2.acom.ucar.edu).  MOSAIC aerosol physics and chemistry are developed by a joint collaboration between PNNL and NCAR and the MOZART gas phase chemistry scheme is developed at NCAR.


%#---- Instructions for running a WRFChem4 Arctic simulation on spirit with
%# MOZART-MOSAIC-AQ-4bin chemistry/aerosols
#
# Louis Marelle, 2023/12/15
#

This runs a low-resolution (100 km) quasi-hemispheric Arctic WRF-Chem
simulation for 2012-02-15 to 2012-02-16, with MOZART-MOSAIC gas-phase+aerosol

# WRF-Chem is the Weather Research and Forecasting model, including chemistry.
# It is a 3D, limited-area (regional) model used for weather forecasts,
# long-term regional climate projections, air-quality forecasts and
# atmospheric process studies.
#
# For more info on WRF and WRF-Chem
# https://www2.mmm.ucar.edu/wrf/users/
# For an introduction to the model,
# https://www.youtube.com/watch?v=wzSu-343b-0
#
# The WRF-Chem model is run in 4 main steps. First, WPS (the WRF preprocessing
# system), then real.exe (program creating the main wrf input and boundary
# files), then additional WRF-Chem preprocessors for emissions and boundary
# conditions, then wrf.exe (the WRF-Chem model). How to run these programs is
# explained below.
#
# Before running this test case, you need to compile the following programs, or a
# compiled version of these programs need to be copied to your own space:
# - WRF-chem model (WRF model with chemistry enabled)
# - WPS (WRF preprocessor) compiled for your WRF version
# - WRF-Chem preprocessors (mozbc, wesely, exo_coldens, megan_bio_emiss, fire_emis)
# You can find compilation instructions in WRF-compile-scripts
#
# I also recommend going through the WRF online tutorial before running this
# test, to better understand how the model works. I also recommend running the
# WRF met-only test in wrf-met/ first (WRF without the chemistry)

Create a directory for storing WRF output on your /data/ space, for example
 mkdir /data/$(whoami)/WRF
 mkdir /data/$(whoami)/WRF/WRF_OUTPUT
