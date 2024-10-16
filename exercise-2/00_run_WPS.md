# First step - running WPS
# *NB: his is the same as we did in exercise-1. No need to run it again if you did exercise-1 already*

WPS is the WRF Pre-Processing System. It comprises several codes that together transform the meteorological input data to a WRF-compatible format.
The scripts and parameter files needed to run WPS for this exercise are located on dahu here: `/bettik/PROJECTS/pr-regionalchem/laperer/WRFChem_training/exercise-1/00_run_WPS/`.
Please copy/paste this folder to your home directory on dahu, from where you will launch the runs.
- `namelist.wps` contains the parameters defining the geographical information on the simulation domain, the projection...
- `jobscript_wps.sh` is the script that will be used to launch the sequence of pre-processors associated with WPS on Dahu
- the other files in the folder are parameter files that do not need to be modified for our exercise

## Input files

Input files from ERA5 and FNL for this exercise have been pre-downloaded so you do not have to deal with the acquisition in this exercise, although this is something you should be able to do in the future. 

> FNL can be downloaded from `https://rda.ucar.edu/datasets/ds083.2/`. 
> Example Python scripts for downloading ERA5 data on sufrace and pressure level needed to drive WPS are provided here as `dl_era5_surface.py` and `dl_era5_pressure_levels.py`. 

On Dahu, the data we have downloaded can be found here: `/bettik/PROJECTS/pr-regionalchem/laperer/WRF_INPUT/met_data/`

Other input files, such as geographical surface information are also already on Dahu at `/bettik/PROJECTS/pr-regionalchem/laperer/WPS_INPUT/geog_wrf4/`

## Paths and namelist options

The `namelist.wps` parameter file defines the domain that will be run including the spatial resolution, the number of grid cells, the nested domains if any, the central lat/lon coordinates, the projection. The path to geographical input data is also given in the file. For the purpose of the exercise there is nothing to modify in this file but please have a look at it.

In the `jobscript_wps.sh` script, you need to update a few paths for where to store the output data (OUTDIR_ROOT, SCRATCHDIR_ROOT, SUBMIT_DIR), and give a name to your simulation case (CASENAME). Typically, you will want OUTDIR_ROOT to point to `/bettik/PROJECTS/pr-regionalchem/username/...`, SCRATCHDIR_ROOT to `/silenus/PROJECTS/pr-regionalchem/username/...` and SUBMIT_DIR should be your current directory (hence the $PWD in the provided script). You can also adjust in this file the run dates and the allocated computing power. For WPS you should run on a single core.


## Launching WPS

When you have setup your `namelist.wps` and adjusted the paths in the `jobscript_wps.sh`, you can launch WPS using the following command: 
> oarsub -S ./jobscript_wps.sh

NB: you might need to `chmod +x jobscript_wps.sh` before launching.

To follow the status of your run, prompt
> oarstat -u username
