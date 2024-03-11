# To watch
We recommend you watch this video before starting this exercise
**Running the WRF Model (for Real and Ideal Cases)**
> https://www.youtube.com/watch?v=yixvMF-g0nc

In this tutorial we will only work with Real cases.

# Exercise obejctive 
Setup of a regional model domain over Europe, run a real case for the dates 
> 1 - 7 January 2019

# WRF initial and boundary conditions 
The WRF model requires initial and boundary conditions to start a simulation. These conditions provide the starting point for the model to simulate the evolution of the atmosphere over time. 

## Initial Conditions
Initial conditions describe the state of the atmosphere at the beginning of the simulation period. They include variables such as temperature, humidity, wind speed, and direction at various levels in the atmosphere. 

## Lateral Boundary Conditions
Boundary conditions specify how the atmosphere interacts with its surroundings at the edges of the simulation domain. They include information about the incoming flow of air, temperature, humidity, and other relevant variables at the boundaries of the model domain. Boundary conditions are often derived from global weather models or from observational data collected at locations surrounding the domain of interest.

## Lower Boundary Conditions
In order to run long simulations with changing surface properties (snow cover, sea ice, etc) we recommend you use the lower boundary condition option available in WRF.  
To adjust or "fix" the lower boundary conditions to be those read in from the analyses/reanalyses in WRF, you typically modify parameters related to land surface properties and processes. This can be done by using WRF lower boundary condition options in the WRF namelist (see below).

## What we use
In this tutorial we will use initial conditions from ERA5 
> https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview
> https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

or from FNL
> https://rda.ucar.edu/datasets/ds083.2/


# Pre-defined model domain
We have pre defined a model domain for you, you can find this on Gricad here

# WRF nudging options 


# Namelists
The main way you will set options and decide how you will use the WRF modeling system and preprocessing system is through namelists.


## What is the WPS namelist?
> 

## What is the WRF namelist?
> The WRF namelist is the input file is used for both the real.exe and wrf.exe executables.  We have provided you an example WRF namelist with the following options chosen:



