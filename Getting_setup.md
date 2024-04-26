The first steps are adapted from [the official GRICAD documentation](https://gricad-doc.univ-grenoble-alpes.fr/en/hpc/connexion/) and [the informal IGE documentation](https://github.com/ige-calcul/public-docs/blob/main/clusters/Gricad/dahu.md). The following ones are tailored for the training purposes. Please refer to the official documentation in case of problems.

1. Getting a [Perseus NG account](https://perseus.univ-grenoble-alpes.fr/) (more details [here](https://gricad-doc.univ-grenoble-alpes.fr/en/services/perseus-ng/account/))

2. Joining the `pr_regionalchem` project on **Perseus**: *Home > My projects > Join a project* (more details [here](https://gricad-doc.univ-grenoble-alpes.fr/en/services/perseus-ng/project/))

3. Use your previously created `<perseus-login>` and its password to check your access to clusters:
    - connect first to the proxy: `ssh <perseus-login>@trinity.u-ga.fr`
    - and then connect to the cluster : `ssh dahu`


4. Set up your SSH key access
    - Generate SSH key pair (or reuse one) in a terminal - you'll be ask to set a password to protect your SSH private key (not needed but recommanded):
      ```
      ssh-keygen
      ```

    - Configure your `.ssh/config` file for SSH easy access to clusters by adding the following lines (**don't forget to change `<perseus-login>` in both place**):
      ```
      Host trinity.u-ga.fr rotule.u-ga.fr access-gricad.u-ga.fr
        User <perseus-login>

      Host dahu.ciment
        Hostname dahu
        User <perseus-login>
        ProxyJump access-gricad
      ```

    - Copy your SSH public key to the cluster proxies and then to the DAHU cluster (password for the SSH private key will be requested first if set, then the Perseus Account password):
      ```
      ssh-copy-id rotule.u-ga.fr
      ssh-copy-id trinity.u-ga.fr
      ssh-copy-id dahu.ciment
      ```

      - Connect to the DAHU cluster using:
      ```
      ssh dahu.ciment
      ```

5. Configure nix environment, python packages and launch a Jupyter notebook on a node (details on nix configuration [here](https://gricad-doc.univ-grenoble-alpes.fr/en/hpc/softenv/nix/))
    - After connecting to dahu compute the following command lines to install Jupyter:
      ```
      source /applis/site/nix.sh
      nix-env --switch-profile $NIX_USER_PROFILE_DIR/jupyter
      ```
      ```
      cat > ~/python.nix <<EOF
      with import <nixpkgs> {};
      pkgs.python310.withPackages (ps: with ps; [ ipython jupyter ])
      EOF
      ```
      ```
      nix-env -i -f ~/python.nix
      ```

    - Set-up python packages:

      In the same terminal than previous commands:
      ```
      nano ~/.config/nixpkgs/config.nix
      ```
      You should open a file structured as follow:
      ```
      # Python environment
      pythonEnv = pkgs.python3.withPackages (ps: with ps; [

            #################################################
            # You can list your python packages below and
            # install a python environement with:
            #      nix-env -f "<nixpkgs>" -iA pythonEnv
            #################################################
            numpy ipython virtualenv pip notebook
      ]
      ```
      Add to the line with package names: "numpy ipython virtualenv pip notebook" the python package you need. For this training, you can add those package:
      ```
      numpy ipython virtualenv pip notebook xarray matplotlib cartopy netcdf4 pyproj shapely pandas
      ```
      Then download your new python environment :
      ```
      nix-env -f "<nixpkgs>" -iA pythonEnv
      ```

    - Jupyter configuration :
      ```
      jupyter notebook --generate-config
      ```
      ```
      echo "c.NotebookApp.open_browser = False" >> ~/.jupyter/jupyter_notebook_config.py
      echo "c.NotebookApp.ip = '0.0.0.0'" >> ~/.jupyter/jupyter_notebook_config.py
      ```
    - Launch a Job on a DAHU node :
      ```
      oarsub -I --project pr-regionalchem -l /core=1,walltime=2:00:00
      ```
      The walltime specify the time of your job. After this time your jupyter notebook will be killed. Change it as you wish.
      ```
      source /applis/site/nix.sh
      jupyter notebook
      ```
      The last command will return usefull information :
      ```
      [I 15:48:01.768 NotebookApp] Writing notebook server cookie secret to /home/user/.local/tmp/jupyter/notebook_cookie_secret
      [I 15:48:09.697 NotebookApp] Serving notebooks from local directory: /home/user
      [I 15:48:09.697 NotebookApp] The Jupyter Notebook is running at:
      [I 15:48:09.697 NotebookApp] http://(dahu38 or 127.0.0.1):8888/?token=524bfbd752f764ef447509bde0a5d8f9169ab603ce33f966
      [I 15:48:09.697 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
      [C 15:48:09.727 NotebookApp]
      To access the notebook, open this file in a browser:
          file:///home/user/.local/tmp/jupyter/nbserver-276827-open.html
      Or copy and paste one of these URLs:
          http://(dahu38 or 127.0.0.1):8888/?token=524bfbd752f764ef447509bde0a5d8f9169ab603ce33f966
      ```
    - Open the Jupyter notebook on your browser. From a local shell :
      ```
      ssh -fNL 8888:node_name:8888 dahu.ciment
      ```
      The name of the node is given in the information of Jupyter notebook command. Here it is dahu38, also change the number 8888 by the one corresponding to the URL.
    - When the ssh tunnel is setup, you can open the Jupyter notebook in your browser with the previous URL with the name of the node replaced by localhost. In this case :
      ```
      http://localhost:8888/?token=524bfbd752f764ef447509bde0a5d8f9169ab603ce33f966
      ```
