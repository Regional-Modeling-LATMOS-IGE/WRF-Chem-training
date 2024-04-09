# Second step - running *real*

*real* is the pre-processor that generates initial and boundary conditions from the forcing data for use in WRF. 

The scripts and parameter files needed to run *real* for this exercise are located on dahu here: `/bettik/PROJECTS/pr-regionalchem/laperer/WRFChem_training/01_run_real/`.
Please copy/paste this folder to your home directory on dahu, from where you will launch the runs.
- `namelist.input` is the parameter file for defining all WRF options. This file should be identical to the one that will be used in the next to run WRF.
- `jobscript_real.sh` is the script that will set the paths, load the libraries and launch `real.exe`
- the other files in the folder are parameter files that do not need to be modified for our exercise

## Input files

Input files ...

## Paths and namelist options

The `namelist.input` parameter file contains all the parameters needed to run the WRF model, including the time controls (simulation dates, restart simulation...), the domain details (number of grid points, resolution... should be consistent with `namelist.wps`), the physics parameterizations to use, parameters for the dynamics and options for handling the boundary conditions.

In the `jobscript_real.sh` script, you need to update a few paths for where to store the output data (OUTDIR_ROOT, SCRATCHDIR_ROOT, SUBMIT_DIR), and include your simulation case (CASENAME - should be the same as in `jobscript_wps.sh`). 
You also need to setup the run dates and the allocated computing power.


## Launching real

When you have setup your `namelist.input` and adjusted the paths in the `jobscript_real.sh`, you can launch real using the following command: 
> oarsub -S ./jobscript_real.sh

To follow the status of your run, prompt
> oarstat -u username
