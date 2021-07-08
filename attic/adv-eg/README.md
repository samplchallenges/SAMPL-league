# Description:
* This container takes a ligand SMILES string and receptor, uses `oechem` to charge the ligand and add 3D coordinates, then uses AutoDock Vina to dock the ligand into the specified binding site. 
* This `README.md` details how to setup, build and run the Autodock Vina docker container. 

# Setup:
1. `mkdir SAMPL-league/attic/adv-base/dependencies`
2. Download Autodock Tools linux x86 .tgz file (`autodock_vina_1_1_2_linux_x86.tgz`) from http://vina.scripps.edu/download.html
3. `mv autodock_vina_1_1_2_linux_x86.tgz SAMPL-league/attic/adv-base/dependencies`
4. `cd SAMPL-league/attic/adv-base/dependencies`
5. `tar -xvf autodock_vina_1_1_2_linux_x86.tgz`
6. `rm autodock_vina_1_1_2_linux_x86.tgz`
7. `mv autodock_vina_1_1_2_linux_x86 adv`
8. Download MGL Tools linux x86 .tar.gz (`mgltools_x86_64Linux2_1.5.6.tar.gz`) from http://mgltools.scripps.edu/downloads
9. `mv mgltools_x86_64Linux2_1.5.6.tar.gz SAMPL-league/attic/adv-base/dependencies`
10. `cd SAMPL-league/attic/adv-base/dependencies`
11. `tar -xvf mgltools_x86_64Linux2_1.5.6.tar.gz`
12. `rm mgltools_x86_64Linux2_1.5.6.tar.gz`
13. `mv mgltools_x86_64Linux2_1.5.6 mgl`
14. `vi mgl/install.sh`
15. Change line 6 `TarDir=` to `TarDir="/opt/app/dependencies/mgl/"`
16. Change line 7 `export MGL_ROOT=""` to `export MGL_ROOT="/opt/app/dependencies/mgl/"`
17. Save install.sh
18. `cp oe_license.txt SAMPL-league/attic/adv-eg`


# Build:
1. `cd SAMPL-league/attic/adv-base`
2. `docker build -t adv-base .`
3. `cd SAMPL-league/attic/adv-wprebuild`
4. `docker build -t adv .`


# Run: 
### Options
```
% docker run -it --rm adv --help
Usage: dock [OPTIONS]

  docks the given smiles string into the receptor within the box specified by
  boxsize and center exhaustiveness does not work at this moment

Options:
  -r, --receptor TEXT           path of receptor PDB to dock the ligand into
  -s, --smiles TEXT             SMILES str of ligand to be docked. quote and
                                add white space at the end "CCC "
  --flex TEXT                   Not working yet. flexible sidechains if any
                                pdb
  --c_x FLOAT                   box center x coordinate; must be used with
                                --c_y/z and --sz_x/y/z
  --c_y FLOAT                   box center y coordinate; must be used with
                                --c_x/z and --sz_x/y/z
  --c_z FLOAT                   box center z coordinate; must be used with
                                --c_x/y and --sz_x/y/z
  --sz_x INTEGER                box size in the x direction; must be used with
                                --c_x/y/z and --sz_y/z
  --sz_y INTEGER                box size in the y direction; must be used with
                                --c_x/y/z and --sz_x/z
  --sz_z INTEGER                box size in the z direction; must be used with
                                --c_x/y/z and --sz_x/y
  --b_xmin FLOAT                box coordinate x min must be used with
                                --b_ymin/zmin/xmax/ymax/zmax
  --b_ymin FLOAT                box coordinate y min must be used with
                                --b_xmin/zmin/xmax/ymax/zmax
  --b_zmin FLOAT                box coordinate z min must be used with
                                --b_xmin/ymin/xmax/ymax/zmax
  --b_xmax FLOAT                box coordinate x max must be used with
                                --b_xmin/ymin/zmin/ymax/zmax
  --b_ymax FLOAT                box coordinate y max must be used with
                                --b_xmin/ymin/zmin/xmax/zmax
  --b_zmax FLOAT                box coordinate z max must be used with
                                --b_xmin/ymin/zmin/xmax/ymax
  -n, --num_modes INTEGER       Number of modes to dock
  -e, --exhaustiveness INTEGER  exhaustiveness of the global search, default=8
  --debug                       prints debug print statements when --debug
                                flag is used
  --output-dir PATH             Output Directory
  --help                        Show this message and exit.
```

### Example run commands
`python ../../ever_given/run.py adv --output-dir <PATH> -s <SMILES_str> -r <receptor_path_from_INPUT_DIR> --sz_x <boxsize_x> --sz_y <boxsize_y> --sz_z <boxsize_z> --c_x <center_x> --c_y <center_y> --c_z <center_z> --output-keys <ligand_key>,<receptor_key>` 
* `python ../../ever_given/run.py adv --file-receptor ../data/adv/4w51-cryo.pdb --smiles "CC(C)Cc1ccccc1 " --sz_x 14 --sz_y 14 --sz_z 14 --c_x -32.355 --c_y 7.263 --c_z 2.207 --output-keys docked_ligand,receptor`
