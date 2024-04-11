# Exercise objective 
Setup of a regional model domain over Europe, run a real `WRF-Chem` case for the dates 1 - 7 January 2019.
For this exercise, all of exercise-1 applies, except we add the `Chem` component of the model. Hereafter are only mentioned the specifics of `Chem`. all information related to WRF is in the exercise-1 folder.

<p align="center">
 <img width="460" src="../exercise-1/First_run_domain.png">
</p>


# WRF-Chem preprocessors

## Initial and boundary conditions 
The WRF-Chem model requires initial and boundary conditions to start a simulation. These conditions provide the starting point for the model to simulate the evolution of the atmosphere over time. 

### Initial Conditions
Initial conditions describe the state of atmospheric composition at the beginning of the simulation period. They include concentrations of gases and aerosolss at various levels in the atmosphere. 

### Lateral Boundary Conditions
Boundary conditions specify how the atmosphere interacts with its surroundings at the edges of the simulation domain. They include information about the incoming flow of gases and aerosols at the boundaries of the model domain. Boundary conditions are often derived from global chemistry-transport models.

### Global runs used for boundary and initial conditions 
In this tutorial we will use initial and boundary conditions from CAM-CHEM 
> https://www.acom.ucar.edu/cesm/subset.shtml

For this exercise, we chose to use the MOZART+MOSAIC+AQCHEM chemistry option. Therefore, these initial and boundary conditions are prepared via the `mozbc` pre-processor.
The `mozbc_mozartmosaic.inp` parameter file provides a mapping between the species provided by CAM-CHEM onto the species defined in WRF-Chem for the chosen chemistry option.


## Anthropogenic emissions
Anthropogenic emissions are taken from the CAMS emissions inventory
> https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-emission-inventories?tab=form

They are prepared with the `.py` script


## Fire emissions
Fire emissions are taken from the Fire INventory from NCAR (FINN), Version 1.5
> https://www.acom.ucar.edu/Data/fire/

They are prepared with the `fire_emis` pre-processor. The parameter file for the mapping depends on the chemistry option choice and is, in our case, `fire_emis_mozartmosaic.inp`.

## Biogenic emissions
Biogenic emissions are taken from the Model of Emissions of Gases and Aerosols from Nature (MEGAN).
> https://www2.acom.ucar.edu/modeling/model-emissions-gases-and-aerosols-nature-megan

These emissions are processed with the `megan_bio_emiss` preprocessor, using the parameter file `megan_bioemiss.inp`

# Namelists
The main way you will set options and decide how you will use the WRF modeling system and preprocessing system is through namelists.

