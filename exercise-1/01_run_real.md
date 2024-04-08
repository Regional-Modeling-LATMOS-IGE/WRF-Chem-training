# Second step - running *real*

*real* is ... 

The scripts and parameter files needed to run *real* for this exercise are located on dahu here: `/bettik/PROJECTS/pr-regionalchem/laperer/WRFChem_training/01_run_real/`.
Please copy/paste this folder to your home directory on dahu, from where you will launch the runs.
- `namelist.input` contains ...
- `jobscript_real.sh` is the script that will ...
- the other files in the folder are parameter files that do not need to be modified for our exercise

## Input files

Input files ...

## Paths

In the `jobscript_real.sh` script, you need to update a few paths for where to store the output data (OUTDIR_ROOT, SCRATCHDIR_ROOT, SUBMIT_DIR), and include your simulation case (CASENAME - should be the same as in `jobscript_wps.sh`). 
You also need to setup the run dates and the allocated computing power.


## Launching real

When you have setup your `namelist.input` and adjusted the paths in the `jobscript_real.sh`, you can launch real using the following command: 
> oarsub -S ./jobscript_real.sh

To follow the status of your run, prompt
> oarstat -u username
