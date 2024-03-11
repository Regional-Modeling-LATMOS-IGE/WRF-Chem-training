# What is WRF?
The Weather Research and Forecasting Model (WRF) is a numerical weather prediction system designed to serve both atmospheric research and operational forecasting/hindcasting applications. WRF-Chem is the Weather Research and Forecasting (WRF) model coupled with Chemistry. WRF-Chem is a 3D, limited-area (regional) model used for weather forecasts,  long-term regional climate projections, air-quality forecasts and
 atmospheric process studies. The model simulates the emission, transport, mixing, and chemical transformation of trace gases and aerosols simultaneously with the meteorology.  

## For more info on WRF and WRF-Chem
>  https://www2.mmm.ucar.edu/wrf/users/
## For an introduction to the model, please watch this before you attend the training
>   https://www.youtube.com/watch?v=wzSu-343b-0

# WRF code distribution
WRF is distrubted at:
> https://github.com/wrf-model/WRF

The current version of WRF and WRF-Chem are distributed together, if you are running WRF vs. WRF-Chem depends on how you compile the model and your WRF namelist (main input file) and other input files provided.
The current of the WRF modeling system - as of 11 March 2024 is 4.5.2.  

All information on the WRF modleing system including user mailing lists, documentation, and support from NCAR the main model developer for WRF dynamics is found at
> https://www2.mmm.ucar.edu/wrf/users/download/get_source.html

WRF-Chem is managed by NOAA and the main page for information is found at
>  https://ruc.noaa.gov/wrf/wrf-chem/ 

## WRF and WRF-Chem are modular modeling systems
There are mutliple options for treating processes in the atmosphere, coupling with the land/ocean/ice surfaces, and atmospheric chemistry (trace gas and aerosols).  We will present one model setup that we frequenly use at IGE/LATMOS.  However, other model setups are possible using other atmopserhic, surface, and atmospheric chemsitry setups.

The specific setup for the chemistry part of WRF-Chem presented here is MOZART-4 gas-phase chemistry
> https://gmd.copernicus.org/articles/3/43/2010/

with MOSAIC 4-bin aerosols, original publication of this aerosol model at

>  https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2007jd008782
The coupling of MOZART-MOSAIC 4 bin aerosol scheme is described at

> *GET DOC WHEN NCAR SITE IS BACK UP*

MOSAIC aerosol physics and chemistry are developed by a joint collaboration between PNNL and NCAR and the MOZART gas phase chemistry scheme is developed at NCAR.


# Steps for running WRF and WRF-Chem

The WRF-Chem model is run in 4 main steps. 
- Running the WRF preprocessor: WPS (the WRF Preprocessing System), which consists of mutiple programs
- Real: Creating the main wrf input and boundary files using real.exe
- **WRF-Chem only:** running all additional WRF-Chem preprocessors for emissions and boundary
 conditions
- Run the WRF model using wrf.exe 

# In order to complete this training
Ensure you have a working / compiled version of the WRF/WRF-Chem model, ensure you have a compiled version of WPS that works for your WRF version.  **WRF-Chem only:** Ensure compiled all WRF-Chem preprocessors (mozbc, wesely, exo_coldens, megan_bio_emiss, fire_emis)

# Recommendations
We recommend using/running the WRF met-only exercises first (WRF without the chemistry) before moving on to running WRF-Chem tests/exercieses.  
