# Introduction to SAMPL Containerized Methods

### Purpose:
In [SAMPL4](https://link.springer.com/article/10.1007%2Fs10822-013-9702-2), we learned that human knowledge can be a key factor influencing the success of a computational drug discovery method. To work around this discovery, we are creating an automated arm of SAMPL challenges to run methods head-to-head without human intervention. To accomplish this, we will containerized [Docker](https://www.docker.com/resources/what-container) methods. 

The following tutorial is meant to teach the basics of building a simple docking pose prediction container using Python code and command line programs (specifically Autodock Vina and MGL Tools). 


### Important Note on "Docker" v. "Docking":
Please note that "Docker" and "docking" are two separate things. 
* **"Docker"** is a program that allows you to containerize methods, essentially taking out human intervention from your containerized program. 
* **"Docking"** describes predicting the structure of a complex, in this case a protein-ligand complex.


### Expected Background Knowledge
* Basic knowledge of Python
* Basic knowledge of Linux/UNIX command line


### Software Requirements
* Linux or Mac operating system
* [Docker Desktop](https://www.docker.com/products/docker-desktop)
* [Docker SDK for Python](https://pypi.org/project/docker/)
* [Python 3](https://www.python.org/downloads/)


### Brief Docker Usage Tutorial
* To build an image, ensure you are in the directory with your Dockerfile and container code, then run `docker build -t <name>:<tag/version> .` (i.e. `docker build -t adv:0.1 .` or `docker build -t adv:latest .`)
* Use the command `docker images` to list out your built images
* To delete Docker images, use `docker images` to list your current images and their IMAGE IDs, then run the command `docker image rm <IMAGE IDs>`


### Pre-Built Autodock Vina Container

A working version of the Autodock Vina container we will build in this tutorial can be found at [Docker hub](https://hub.docker.com/repository/docker/osatom/adv-tutorial). To play with this container, please use the following steps: 
1. Pull the "adv-tutorial" docker container: `docker pull osatom/adv-tutorial:latest`
2. Change directories into the "examples" directory
3. Run the command: `python ever_given/run.py osatom/adv-tutorial:latest --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`
4. The results will be stored in the directory "examples/evergiven_output"


# Tutorial: Build an AutoDock Vina Containerized Method

### Outline: 
* Section 1: Build the Autodock Vina base container
  * 1.1 Setup
  * 1.2 Create a conda environment
  * 1.3 Create a Dockerfile
  * 1.4 Add the Autodock Vina and MGL Tools executables
  * 1.5 Update the Dockerfile to include Autodock Vina and MGL Tools installations
  * 1.6 Build the base container
* Section 2: Build the Autodock Vina Docing methods container
  * 2.1 Setup


### Section 1: Build the Autodock Vina base container
In this section, we will base container that has all necessary packages and programs installed that will not be or rarely be changed to improve build time. This way, the container we write our docking code in can be quickly built as we develop it. 


**1.1: Setup**
1. Create a directory called "adv-tutorial-base": `mkdir adv-tutorial-base`
2. Change directories to "adv-tutorial-base": `cd adv-tutorial-base`


**1.2: Create a conda environment**

In 1.2, we will run the "continuumio/miniconda3" container to dynamically create the conda environment we need.
1. Open a terminal
2. Start the container: `docker run -it --rm continuumio/miniconda3` upon running this command, your command line prompt should change
3. Create a conda env called advenv: `conda create --name advenv`
4. Activate advenv: `conda activate advenv`
5. Install rdkit: `conda install -c conda-forge rdkit`
6. Install mdtraj: `conda install -c conda-forge mdtraj`
7. Install click: `pip install click`
8. Export the environment: `conda env export -n advenv`
9. Copy the output from the export command in step 7 to be pasted into a file in the next section.
10. Exit the container: `exit`
11. Create and open a file called "environment.yml" and paste the output you previously copied at the end of 1.2
12. Change the first line of the file "name: advenv" to "name: base"
13. Delete the last line of the file: "prefix: /opt/conda/envs/advenv"
14. Save the changes to "environment.yml" and exit


**1.3: Create a Dockerfile**

In 1.3, we will begin creating a Dockerfile which contains the instructions required to build the container.
1. Create and open a file called "Dockerfile"
2. Copy the following lines into Dockerfile
  ```
  FROM continuumio/miniconda3:4.9.2-alpine  
  # tells the container to inherit from a miniconda container

  WORKDIR /opt/app/   # set the work directory

  COPY . ./    #  copy all the files and directories into the container

  RUN conda env update -f environment.yml && \
      conda clean --all --yes      # install the packages in environment.yml into containers

  ENV PATH="/root/.local/bin:$PATH"      # set the path
```
3. Save the changes to Dockerfile and exit


**1.4: Add the Autodock Vina and MGL Tools executables**

In 1.4, we will incorporate Autodock Vina and MGL Tools into our base container. Please do not change which installer you download based on your native operating system. These installers will be used inside the docker container which uses a Linux x86 system.
1. Create a directory called dependencies: `mkdir dependencies`
2. Download Autodock Tools linux x86 "autodock_vina_1_1_2_linux_x86.tgz" from [http://vina.scripps.edu/download.html]
3. Untar "autodock_vina_1_1_2_linux_x86.tgz": `tar -xvf dependencies/autodock_vina_1_1_2_linux_x86.tgz`
4. Delete the .tgz file: `rm autodock_vina_1_1_2_linux_x86.tgz`
5. Rename "autodock_vina_1_1_2_linux_x86/" to "adv/": `mv autodock_vina_1_1_2_linux_x86 adv`
6. Move the "adv" directory inside the "dependencies" directory: `mv adv dependencies`
7. Download MGL Tools linux x86 `mgltools_x86_64Linux2_1.5.6.tar.gz` from [http://mgltools.scripps.edu/downloads]
8. Untar "mgltools_x86_64Linux2_1.5.6.tar.gz": `tar -xvf mgltools_x86_64Linux2_1.5.6.tar.gz`
9. Delete the .tgz file: `rm mgltools_x86_64Linux2_1.5.6.tar.gz`
10. Rename "mgltools_x86_64Linux2_1.5.6.tar.gz" to "mgl": `mv mgltools_x86_64Linux2_1.5.6.tar.gz mgl`
11. Move "mgl" directory inside the "dependencies" directory: `mv mgl dependencies/`
12. Open "dependencies/mgl/install.sh"
13. Change line 6 from `TarDir=""` to `TarDir="/opt/app/dependencies/mgl/"`
14. Change line 7 from `export MGL_ROOT=""` to `export MGL_ROOT="/opt/app/dependencies/mgl/"`
15. Close and save "dependencies/mgl/install.sh"

**1.5 Update the Dockerfile to include Autodock Vina and MGL Tools installations**
1. Open the Dockerfile and paste the following lines 
  ```
  RUN /opt/app/dependencies/mgl/install.sh

  ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"

  RUN /opt/app/dependencies/adv/bin/vina --help
  ```
2. Save the changes to the `Dockerfile` and exit

**1.6: Build the base container**
1. Build the base container: `docker build -t adv-tutorial-base .`


# Section 2: Build the container with Autodock Vina Docing methods

### Part 1: Write your methods
1. Your method you run as your main should include the following flags to deal with our command line inputs. We typically use `@click.option` to handle this in python
    * `--receptor`: a pdb of the receptor to dock into
    * `--smiles`: SMILES string of the ligand to dock
    * `--hint`: a pdb of a receptor and ligand complex
    * `--hint_radius`: radius from the ligand considered part of the binding site
    * `--hint_molinfo`: ligand resname used in the hint pdb
    * `--output-dir`: directory to output the output files to
2. Ensure your method writes out to the following files:
    * A receptor pdb in the same frame of reference as your docked ligand
    * A ligand pdb, mol2, or sdf in the same frame of reference as your output pdb
3. Ensure your method prints the following lines with `docked_ligand` and `receptor` both no caps exactly as is. This is extremely important as these keys will be used to parse your output to find the file information
    * `docked_ligand {path_to_ligand_file}`
    * `receptor {path_to_receptor_file}`
4. If you would like to decline to predict a pose for a ligand please print the following. This tells us your program is intentionally avoiding these predictions
    * `docked_ligand no_prediction`
    * `receptor no_prediction`

### Part 2: Create your setup.py file
1. Create a file called `setup.py` with the following
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
2. In `py_modules` add any modules you have built that your program requires
3. In `install_requires` add any pip installable packages you did not add in the previous base build that you’d like to add here. You can also go back and rebuild the base build with the new package edit
4. Add an `entry_point`: `{command-to-call-in-Dockerfile}={py_module_with_main}:{function_to_run}`

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
* To use `ever_given/run.py`:
    * any CLI options that are files need to be pre-pended by `--file-`
        * `--file-receptor`
        * `--file-hint`
        * Any file options for ever_given besides `--output-dir` should use underscores rather than dashes `-`

    * You will need to specify the expected output keys as well. For docking it is as follows
        * `--output-keys docked_ligand,receptor`
* `ever_given` will handle all the overhead of running the docker container and mounting the directories so files can be passed between your computer and container
* To run this container from the `examples` directory: `python ever_given/run.py adv --file-receptor <receptor_pdb> --file-hint <hint_pdb> --hint_radius <float> --hint_molinfo <str> --output-keys docked_ligand,receptor`
  * `python ever_given/run.py adv --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`





### Helpful tidbits
* To upload your own container to docker.io

* Eventually, upon building enough Docker images, you may begin to run out of memory. Please remember to regularly delete any Docker images you no longer need. 

* For more detailed tutorials on how to use docker please see the following resources:
  * Official Docker documentation: https://docs.docker.com/get-started/
  * Brief Docker Tutorial (12m): https://www.youtube.com/watch?v=YFl2mCHdv24
  * Docker Beginners Course (2hrs): https://www.youtube.com/watch?v=fqMOX6JJhGo
  * Brief video on inner workings of Docker (15m): https://www.youtube.com/watch?v=rOTqprHv1YE 
  * YouTube is a great resource for learning Docker, feel free to search for other tutorials that suit your specific needs as well 
