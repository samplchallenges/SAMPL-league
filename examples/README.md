# Introduction to SAMPL Containerized Methods

## Directory Structure
  ```
examples
├── CondaEnvInstructions.pdf
├── ContainerRequirements.md
├── README.md
├── adv
│   ├── Dockerfile
│   ├── autodock.py
│   ├── run_autodock.py
│   └── setup.py
├── adv-base
│   ├── Dockerfile
│   ├── README.md
│   ├── environment.yml
│   └── setup.py
├── data
│   ├── hint.pdb
│   └── receptor.pdb
└── ever_given
    ├── LICENSE
    ├── Pipfile
    ├── Pipfile.lock
    ├── README.md
    ├── ever_given
    │   └── wrapper.py
    ├── ever_given.egg-info
    ├── pyproject.toml
    └── run.py
    ├── setup.cfg
    ├── setup.py
    └── tests
  ```

## Tutorial Purpose:
The following tutorial is meant to teach you the basics of building a basic docking pose prediction container using python code and command line programs (specifically Autodock Vina and MGL Tools). 

For a more generalized explanation about inputs and outputs please see `ContainerRequirements.md`

For a more generalized explanation on how to create a conda environment inside your container please see `CondaEnvInstructions.pdf`


## Important Note on Docker v. Docking:
Please note that Docker and docking are two separate things. 
* **Docker** is a program that allows you to containerize methods, essentially taking out human intervention from your containerized program. 
* **Docking** describes predicting the structure of a complex, in this case a protein-ligand complex.

## Tutorial Requirements
* Download Docker Desktop: https://www.docker.com/products/docker-desktop
* Open the Docker Desktop app to start the Docker daemon. 
  * If you ever have trouble running Docker commands from the command line, please make sure you have already started your Docker Desktop app
* From the command line, run `pip install docker` 
* Go to Docker hub https://hub.docker.com/ and create an account to store your Docker containers
  * Please remember a free account can only have one private container at a time, so if your container has sensitive information be sure you are uploading it to your one private container


## Brief Docker Usage Tutorial
* Context: Every time you build your program/Docker container a "Docker image" of your program is created. This image is built by using the build instructions in the `Dockerfile` (explained further below)
* To build an image, ensure you are in the directory with your `Dockerfile` and container code. Then run `docker build -t <name>:<tag/version> .` (i.e. `docker build -t adv:0.1 .` or `docker build -t adv:latest .`
* From the command line, run: `docker login` and follow the prompts to log in to the Docker hub account you created in the previous section. This will allow you to 
push docker images directly to your Docker hub account from the command line
* To push a container to your Docker hub use the command: `docker push <name>:<tag/version>`
* Use the command `docker images` to list out your built images
* Eventually, your Docker images will build up and you may run out of memory. To delete docker images, use `docker images` to list your current images and their IMAGE IDs, then run the command `docker image rm <IMAGE IDs>`
* For more detailed tutorials on how to use docker please see the following resources:
  * Official Docker documentation: https://docs.docker.com/get-started/
  * Brief Docker Tutorial (12m): https://www.youtube.com/watch?v=YFl2mCHdv24
  * Docker Beginners Course (2hrs): https://www.youtube.com/watch?v=fqMOX6JJhGo
  * Brief video on inner workings of Docker (15m): https://www.youtube.com/watch?v=rOTqprHv1YE 


# Tutorial: Build an AutoDock Vina Containerized Method

## Build your base container
We want to build a base container that has all necessary packages and programs installed that will not be or rarely be changed. This way, the container you write your code in can inherit from this pre-built container, and your container will build much faster.

For a more generalized explanation on how to create a conda environment inside your container please see `CondaEnvInstructions.pdf`

### Part 1: Create the conda environment
In this section, we will run the `continuumio/miniconda3` container to dynamically create the conda environment we need.
1. Open a terminal
2. Start the container: `docker run -it --rm continuumio/miniconda3` upon running this command, your command line prompt should change
3. Create a conda env called advenv: `conda create --name advenv`
4. Activate advenv: `conda activate advenv`
5. Install rdkit: `conda install -c conda-forge rdkit`
6. Install mdtraj: `conda install -c conda-forge mdtraj`
7. Install click: `pip install click`
8. Export the environment: `conda env export -n advenv`
9. Copy the output from the export command in step 7
10. Exit the container: `exit`

### Part 2: Create a Dockerfile with build instructions and add the lines to install the conda environment
1. Navigate to the directory you will begin writing your container
2. Create and open a file called `environment.yml` and paste the output you previously copied
3. Change the first line of the file `name: advenv` to `name: base`
4. Delete the last line: `prefix: /opt/conda/envs/my-rdkit-env`
5. Create and open a file called `Dockerfile` which will contain your build instructions
6. Copy the following lines into the Dockerfile:
  ```
  FROM continuumio/miniconda3:4.9.2-alpine  
  # tells the container to inherit from a miniconda container

  WORKDIR /opt/app/   # set the work directory

  COPY . ./    #  copy all the files and directories into the container

  RUN conda env update -f environment.yml && \
      conda clean --all --yes      # install the packages in environment.yml into containers

  ENV PATH="/root/.local/bin:$PATH"      # set the path
```

### Part 3: Add the Autodock Vina executables and MGL directories
In this section, please do not change which installer you download based on your native operating system. These installers will be used inside the docker container which uses a Linux x86 system.
1. Create a directory called dependencies: `mkdir dependencies`
2. Download Autodock Tools linux x86 .tgz file (`autodock_vina_1_1_2_linux_x86.tgz`) from http://vina.scripps.edu/download.html
3. Move `autodock_vina_1_1_2_linux_x86.tgz` to the `dependencies` directory, untar and rename the directory to `adv`; remove the .tgz file
4. Download MGL Tools linux x86 .tar.gz (`mgltools_x86_64Linux2_1.5.6.tar.gz`) from http://mgltools.scripps.edu/downloads
5. Move `mgltools_x86_64Linux2_1.5.6.tar.gz` to the `dependencies` directory, untar and rename the directory to `mgl`; remove the tar file
6. Open `mgl/install.sh`
7. Change line 6 `TarDir=` to `TarDir="/opt/app/dependencies/mgl/"`
8. Change line 7 export MGL_ROOT="" to export MGL_ROOT="/opt/app/dependencies/mgl/"
9. Close and save `mgl/install.sh`
10. Add the following to your `Dockerfile`
  ```
  RUN /opt/app/dependencies/mgl/install.sh

  ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"

  RUN /opt/app/dependencies/adv/bin/vina --help
  ```

### Part 4: Build your base container
1. Run the following command: docker build -t adv-base .


## Build your container with your methods

### Part 1: Write your methods
1. Your method you run as your main should include the following flags to deal with our command line inputs. We typically use `@click.option` to handle this in python
    * `--receptor`: a pdb of the receptor to dock into
    * `--smiles`: SMILES strof the ligand to dock
    * `--hint`: a pdb of a receptor and ligand complex
    * `--hint_radius`: radius from the ligand considered part of the binding site
    * `--hint_molinfo`: ligand resname used in the hint pdb
    * `--output-dir`: directory to output the output files to
2. Ensure your method writes out to the following files:
    * A receptor pdb in the same frame of reference as your docked ligand
    * A ligand pdb, mol2, or sdf in the same frame of reference as your output pdb
3. Ensure your method prints the following lines with `docked_ligand` and `receptor` both no caps exactly as is. This is extremely important as these are the keys we will look for in your output to find the file information
    * `docked_ligand {path_to_ligand_file}`
    * `receptor {path_to_receptor_file}`
4. If you would like to decline to predict a pose for a ligand please print the following. This tells us your program is intentionally avoiding these predictions
    * `docked_ligand no_prediction`
    * `receptor no_prediction`

### Part 2: Create your setup.py file
1. Create a file called setup.py with the following
  ```
  from setuptools import setup

  setup(
      name='AutoDock-rdkit',
      version='0.1',
      py_modules=[
  		    ‘run_autodock’,
  	  ],
      install_requires=[
          'Click',
      ],
      entry_points='''
          [console_scripts]
          run-autodock=run_autodock:autodock
      '''
  )
  ```
2. In py_modules add any modules you have built that your program requires
3. In install_requires add any pip installable packages you did not add in the previous base build that you’d like to add here. You can also go back and rebuild the base build with the new package edit
4. Add an entry_point: {command-to-call-in-Dockerfile}={py_module_with_main}:{function_to_run}

### Part 3: Create your Dockerfile
1. Create a file called Dockerfile and copy the following
  ```
  FROM adv-base

  WORKDIR /opt/app/

  COPY setup.py ./

  RUN pip install .

  COPY run_autodock.py ./

  ENV PATH="/root/.local/bin:$PATH"

  RUN ls -l /opt/app

  ENTRYPOINT ["run-autodock"]
  ```

2. In the `FROM` inherit from your previous base build (see section "Build your base container")
3. Set the `WORKDIR`
4. Copy your `setup.py` and other py modules into the container
5. Install `setup.py` using `RUN pip install .`
6. Set the `ENV PATH`
7. Add the `ENTRYPOINT` you declared in step 4 of the previous section

### Part 4: Build your container
1. Build the container: `docker build -t adv .`


## Test/Run your container
To run your container, and ensure that it will work with our infrastructure, please use `ever_given`. `ever_given` handles file input/output and volume mounting
* In the `ever_given` directory, there is file called `run.py`, this is the executable that mimics our infrastructure to run your container
* To use ever_given/run.py:
    * any CLI options that are files need to be pre-pended by `--file-`
        * `--file-receptor`
        * `--file-hint`
        * Any file options for ever_given besides `--output-dir` should use underscores rather than dashes `-`

    * You will need to specify the expected output keys as well. For docking it is as follows
        * `--output-keys docked_ligand,receptor`
* `ever_given` will handle all the overhead of running the docker container and mounting the directories so files can be passed between your computer and container
* To run this container from the `examples` directory: `python ever_given/run.py adv --file-receptor <receptor_pdb> --file-hint <hint_pdb> --hint_radius <float> --hint_molinfo <str> --output-keys docked_ligand,receptor`
  * `python ever_given/run.py adv --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`
