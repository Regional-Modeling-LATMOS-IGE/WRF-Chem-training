# Running WRF-Chem, finally!

The scripts and parameter files needed to run *WRF-Chem* for this exercise are located on dahu here: `/bettik/PROJECTS/pr-regionalchem/laperer/WRFChem_training/exercise-2/02_run_WRFChem/`.
Please copy/paste this folder to your home directory on dahu, from where you will launch the runs.
- `namelist.input.2019` is the parameter file for defining all WRF-Chem options.
- `jobscript_wrf.sh` is the script that will set the paths, load the libraries and launch `wrf.exe`

## Paths and namelist options

The `namelist.input.2019` parameter file contains all the parameters needed to run the *WRF-Chem* model, including the time controls (simulation dates, restart simulation...), the domain details (number of grid points, resolution... should be consistent with `namelist.wps`), the physics parameterizations to use, parameters for the dynamics, options for handling the boundary conditions, and chemistry options.

In the `jobscript_wrf.sh` script, you need to update a few paths for where to store the output data (OUTDIR_ROOT, SCRATCHDIR_ROOT, SUBMIT_DIR), and include your simulation case (CASENAME - should be the same as in `jobscript_wps.sh` and `jobscript_real.sh`). You also need to setup the run dates and the allocated computing power.


## Launching WRF-Chem

When you have setup your `namelist.input.2019` and adjusted the paths in the `jobscript_wrf.sh`, you can launch *WRF-Chem* using the following command: 
> oarsub -S ./jobscript_wrf.sh

NB: you might need to `chmod +x jobscript_wrf.sh` before launching.

To follow the status of your run, prompt
> oarstat -u username
