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
In the context of the WRF model, the 'wrf_low_input' file contains specifications for the lower boundary conditions, particularly for the land surface processes. To adjust or "fix" the lower boundary conditions in WRF, you typically modify parameters related to land surface properties and processes. Here's how you can do it:

Edit the wrf_low_input file: This file contains settings related to land surface properties such as land use, vegetation type, soil type, and land cover. You can adjust these parameters based on the characteristics of the region you are simulating.

Land Use/Land Cover (LU/LC) Categories: Ensure that the land use and land cover categories in the wrf_low_input file accurately represent the area being simulated. These categories influence surface energy and moisture fluxes, which in turn affect the lower boundary conditions.

Soil Parameters: Adjust soil parameters such as soil texture, moisture content, and thermal properties to reflect the characteristics of the soil in the region of interest. These parameters affect soil temperature, moisture, and heat fluxes at the land surface.

Vegetation Parameters: If your simulation includes vegetation, ensure that vegetation parameters such as leaf area index (LAI), roughness length, and stomatal resistance are appropriate for the vegetation type and seasonal variations in the area being simulated.

Urban Parameters (if applicable): If your simulation includes urban areas, make sure to specify urban parameters such as building height, albedo, and roughness length to represent the urban surface accurately.

Surface Fluxes: Adjust parameters related to surface fluxes such as roughness length, albedo, and emissivity to account for the effects of surface properties on heat, moisture, and momentum exchange between the land surface and the atmosphere.

Sensitivity Analysis: Conduct sensitivity tests by varying key parameters in the wrf_low_input file to understand their impact on model results. This can help you identify which parameters have the most significant influence on the lower boundary conditions and how they affect the simulation.

After making adjustments to the wrf_low_input file, rerun the WRF model to see how the changes impact the simulation results. It's essential to validate the model output against observational data to ensure that the adjusted lower boundary conditions produce realistic simulations.

## What we use
In this tutorial we will use initial conditions from ERA5 
> https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview
> https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

or from FNL
> https://rda.ucar.edu/datasets/ds083.2/


# Pre defined model domain
We have pre defined a model domain for you, you can find this on Gricad here

# WRF nudging options 


# Namelists
The main way you will set options and decide how you will use the WRF modeling system and preprocessing system is through namelists.


## What is the WPS namelist?
> 

## What is the WRF namelist?
> The WRF namelist is the input file is used for both the real.exe and wrf.exe executables.  We have provided you an example WRF namelist with the following options chosen:



