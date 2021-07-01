# Setup:
1. `mkdir SAMPL-league/attic/adv-wprebuild/dependencies`
2. Download Autodock Tools linux x86 .tgz file (`autodock_vina_1_1_2_linux_x86.tgz`) from http://vina.scripps.edu/download.html
3. `mv autodock_vina_1_1_2_linux_x86.tgz SAMPL-league/attic/adv-wprebuild/dependencies`
4. `cd SAMPL-league/attic/adv-wprebuild/dependencies`
5. `tar -xvf autodock_vina_1_1_2_linux_x86.tgz`
6. `rm autodock_vina_1_1_2_linux_x86.tgz`
7. `mv autodock_vina_1_1_2_linux_x86 adv`
8. Download MGL Tools linux x86 .tar.gz (`mgltools_x86_64Linux2_1.5.6.tar.gz`) from http://mgltools.scripps.edu/downloads
9. `mv mgltools_x86_64Linux2_1.5.6.tar.gz SAMPL-league/attic/adv-wprebuild/dependencies`
10. `cd SAMPL-league/attic/adv-wprebuild/dependencies`
11. `tar -xvf mgltools_x86_64Linux2_1.5.6.tar.gz`
12. `rm mgltools_x86_64Linux2_1.5.6.tar.gz`
13. `mv mgltools_x86_64Linux2_1.5.6 mgl`
14. `vi mgl/install.sh`
15. Change line 6 `TarDir=` to `TarDir="/opt/app/dependencies/mgl/"`
16. Change line 7 `export MGL_ROOT=""` to `export MGL_ROOT="/opt/app/dependencies/mgl/"`
17. Save install.sh
18. `cp oe_license.txt SAMPL-league/attic/adv-wprebuild`


# Build:
1. `cd SAMPL-league/attic/adv-base`
2. `docker build -t adv-base:0.1 .`
3. `cd SAMPL-league/attic/adv-wprebuild`
4. `docker build -t adv .`


# Run: 
### Options
```
docker run -it --rm adv-pb --help
Usage: dock [OPTIONS]

  docks the given smiles string into the receptor within the box specified by
  boxsize and center exhaustiveness does not work at this moment

Options:
  -r, --receptor TEXT             path of receptor PDB to dock the ligand into
  -s, --smiles TEXT               SMILES str of ligand to be docked. quote and
                                  add white space at the end "CCC "
  --boxsize <INTEGER INTEGER INTEGER>...
                                  The size of the box to dock into. Must be
                                  used with --center and without --boxcoords
  --center <FLOAT FLOAT FLOAT>...
                                  The center of the box to dock into. Must be
                                  used with --boxsize and without --boxcoords
  --boxcoords <FLOAT FLOAT FLOAT FLOAT FLOAT FLOAT>...
                                  The minimum and maximum corners of the box
                                  to dock into. Must be used without --boxsize
                                  and --center
  -e, --exhaustiveness TEXT       Not working yet
  --bind_out TEXT                 Directory in the container the outputs are
                                  bound to
  --bind_in TEXT                  Directory in the container the inputs are
                                  bound to
  --help                          Show this message and exit.
  ```

### Example run commands
`docker run -it --rm -v <INPUT_DIR>:<BIND_IN> -v <OUTPUT_DIR>:<BIND_OUT> adv --bind_in <BIND_IN> --bind_out <BIND_OUT> -s <SMILES_str> -r <receptor_path_from_INPUT_DIR> --boxsize <boxsize_x boxsize_y boxsize_z> --center <center_x center_y center_z> ` 


`docker run -it --rm -v <INPUT_DIR>:<BIND_IN> -v <OUTPUT_DIR>:<BIND_OUT> adv --bind_in <BIND_IN> --bind_out <BIND_OUT> -s <SMILES_str> -r <receptor_path_from_INPUT_DIR> --boxcoords <xmin ymin zmin xmax ymax zmax>`
