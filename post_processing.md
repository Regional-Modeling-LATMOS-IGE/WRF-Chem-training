# Post processing WRF outputs on Dahu

- `source /applis/environments/conda.sh`
- `conda create -n myenv python=3.9`
- `conda activate myenv`
- `conda install conda-forge::ncview`

- set up conda environment / jupyter notebook (https://gricad-doc.univ-grenoble-alpes.fr/notebook/hpc/)

Note: to allow ncview to pop up, you need to connect using "ssh -A -Y username@trinity.univ-grenoble-alpes.fr"

# getting an interactive node for plotting
``` oarsub -I --project pr-regionalchem -l /core=1,walltime=2:00:00 ```
