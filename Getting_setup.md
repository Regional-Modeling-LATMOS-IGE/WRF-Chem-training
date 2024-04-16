1. Getting a gricad account: https://gricad-doc.univ-grenoble-alpes.fr/en/services/perseus-ng/account/
  
2. Joining the pr_regionalchel project: https://gricad-doc.univ-grenoble-alpes.fr/en/services/perseus-ng/project/

3. Follow the instructions for getting setup:
- https://gricad-doc.univ-grenoble-alpes.fr/en/hpc/connexion/

4. Checking access:
   Once you created you perseus account, you should have your "perseus-login" which will be used for configuring your access in the next steps
- Make sure that you have access to the cluster
   - connect to the “bastion” : ssh  perseus-login@trinity.u-ga.fr
   - and then connect to a the frontnode of the cluster : ssh f-dahu
- For more information on connecting and doing scientific calculations on dahu see https://github.com/IGE-numerique/ige-calcul/blob/main/clusters/Gricad/dahu.md 


5. Configure your SSH key access setup (Documentation : https://gricad-doc.univ-grenoble-alpes.fr/en/hpc/connexion/)

  - Generate ssh key pair (You only need to do it once) :
      ```
      ssh-keygen -t rsa
      ```
      Don't put any passwords, just tape the Enter command twice
    
  - Copy your public key to rotule and trinity :
      ```
      ssh-copy-id  perseus-login@rotule.u-ga.fr
      ```
      ```
      ssh-copy-id  perseus-login@trinity.u-ga.fr
      ```
  - Configure your config file for ssh connection to dahu :

      Add the following lines to your .ssh/config file
      ```
      Host *
        ServerAliveInterval 30

      Host *.ciment
        User perseus-login
        ProxyCommand ssh -q perseus-login@access-gricad.univ-grenoble-alpes.fr "nc -w 60 `basename %h .ciment` %p"
      ```
  - Then you can add your public key to dahu :
    ```
    ssh-copy-id  dahu.ciment
    ```
  - You should now be able to connect to dahu using :
    ```
    ssh dahu.ciment
    ```

6. Launch a jupyter notebook on a node

  - Jupyter installation : After connecting to dahu compute the following command lines :
    ```
    source /applis/site/nix.sh
    ```
    ```
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
  - Set-up python packages (see [documentation](https://gricad-doc.univ-grenoble-alpes.fr/en/hpc/softenv/nix/)):
    In the same terminal than previous commands,
    ```
    nano ~/.config/nixpkgs/config.nix
    ```
    You should open a file structured as follow :
    ```
    # Python environment
    pythonEnv = pkgs.python3.withPackages (ps: with ps; [

          #################################################
          # You can list your python packages below and
          # install a python environement with:
          #      nix-env -f "<nixpkgs>" -iA pythonEnv
          #################################################
          numpy ipython virtualenv pip notebook
    ```
    Add to the line with package names : "numpy ipython virtualenv pip notebook" the python package you need. For this training, you can add those package :
    ```
    numpy ipython virtualenv pip notebook xarray matplotlib cartopy
    ```
    Then download your new python environement :
    ```
    ix-env -f "<nixpkgs>" -iA pythonEnv
    ```
  - Jupyter configuration :
    ```
    jupyter notebook --generate-config
    ```
    ```
    echo "c.NotebookApp.open_browser = False" >> ~/.jupyter/jupyter_notebook_config.py
    ```
    ```
    echo "c.NotebookApp.ip = '0.0.0.0'" >> ~/.jupyter/jupyter_notebook_config.py
    ```
- Launch a Job on a Gricad node :
    ```
    oarsub -I --project pr-regionalchem -l /core=1,walltime=2:00:00
    ```
    The walltime specify the time of your job. After this time your jupyter notebook will be killed. Change it as you wish.
    ```
    source /applis/site/nix.sh
    ```
    ```
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
- Open the jupyter notebook on your browser. From a local shell :
    ```
    ssh -fNL 8888:node_name:8888 dahu.ciment
    ```
    The name of the node is given in the information of jupyter notebook command. Here it is dahu38, also change the number 8888 by the one corresponding to the URL.
- When the ssh tunnel is setup, you can open the jupyter notebook in your browser with the previous URL with the name of the node replaced by localhost. In this case :
    ```
    http://localhost:8888/?token=524bfbd752f764ef447509bde0a5d8f9169ab603ce33f966
    ```
      
