# First step - running WPS

WPS is the WRF Pre-Processing System. It comprises several codes that together transform the meteorological input data to a WRF-compatible format.
The scripts and parameter files needed to run WPS for this exercise are located on dahu here: `/bettik/PROJECTS/pr-regionalchem/laperer/WRFChem_training/00_run_WPS/`.
Please copy/paste this folder to your home directory on dahu, from where you will launch the runs.
- `namelist.wps` contains the parameters defining the geographical information on the simulation domain, the projection...
- `jobscript_wps.sh` is the script that will be used to launch the sequence of pre-processors associated with WPS on Dahu
- the other files in the folder are parameter files that do not need to be modified for our exercise

## Input files

Input files from ERA5 and FNL for this exercise have been pre-downloaded so you do not have to deal with the acquisition in this exercise, although this is something you should be able to do in the future. FNL can be downloaded from `https://rda.ucar.edu/datasets/ds083.2/`. Example Python scripts for downloading ERA5 data on sufrace and pressure level needed to drive WPS are provided here as `dl_era5_surface.py` and `dl_era5_pressure_levels.py`. 

The data can be found here: `/bettik/PROJECTS/pr-regionalchem/laperer/WRF_INPUT/met_data/`

Other input files, such as geographical surface information are also already on Dahu at `/bettik/PROJECTS/pr-regionalchem/laperer/WPS_INPUT/geog_wrf4/`

## Paths

In the `jobscript_wps.sh` script, you need to update a few paths for where to store the output data (OUTDIR_ROOT, SCRATCHDIR_ROOT, SUBMIT_DIR), and give a name to your simulation case (CASENAME). You can also adjust in this file the run dates and the allocated computing power.


## Launching WPS

When you have setup your `namelist.wps` and adjusted the paths in the `jobscript_wps.sh`, you can launch WPS using the following command: 
> oarsub -S ./jobscript_wps.sh

To follow the status of your run, prompt
> oarstat -u username
