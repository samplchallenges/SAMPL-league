# Setup
1. `cp oe_license.txt SAMPL-league/attic/oedock`

# Build
1. `cd SAMPL-league/attic/oedock-base`
2. `docker build -t oedock-base .`
3. `cd SAMPL-league/attic/oedock`
4. `docker build -t oedock`

# Run

### Options
```
docker run -it --rm oedock --help
Usage: oedock [OPTIONS]

Options:
  -s, --smiles TEXT               SMILES string  [required]
  -r, --receptor TEXT             receptor PDB  [required]
  --boxcoords <FLOAT FLOAT FLOAT FLOAT FLOAT FLOAT>...
                                  The minimum and maximum values of the
                                  coordinates of the box representing the
                                  binding site. Enter in the order: xmin ymin
                                  zmin xmax ymax zmax
  --boxsize <INTEGER INTEGER INTEGER>...
                                  Size of the box to dock into. Enter in the
                                  order of: size_x size_y size_z
  --center <FLOAT FLOAT FLOAT>...
                                  Center of the box to dock into. Enter in the
                                  order of: center_x center_y center_z
  --bind_out TEXT                 directory in the container the outputs are
                                  bound to  [required]
  --bind_in TEXT                  directory in the container the inputs are
                                  bound to  [required]
  --help                          Show this message and exit.
```



### Example run commands
`docker run -it --rm -v <INPUT_DIR>:<BIND_IN> -v <OUTPUT_DIR>:<BIND_OUT> oedock --bind_in <BIND_IN> --bind_out <BIND_OUT> -s <SMILES_str> -r <receptor_path_from_INPUT_DIR> -boxsize <boxsize_x boxsize_y boxsize_z> -center <center_x center_y center_z>`
* Ex: `docker run -it --rm -v $(pwd)/receptors:/data/in -v $(pwd):/data/out oedock --bind_in /data/in --bind_out /data/out -s "CCCc1ccccc1" -r 4w51-cryo.pdb --boxsize 14 14 14 --center -32.355 7.263 2.207`



`docker run -it --rm -v <INPUT_DIR>:<BIND_IN> -v <OUTPUT_DIR>:<BIND_OUT> oedock --bind_in <BIND_IN> --bind_out <BIND_OUT> -s <SMILES_str> -r <receptor_path_from_INPUT_DIR> --boxcoords <xmin ymin zmin xmax ymax zmax>`
* Ex: `docker run -it --rm -v $(pwd)/receptors:/data/in -v $(pwd):/data/out oedock --bind_in /data/in --bind_out /data/out -s "CCCc1ccccc1" -r 4w51-cryo.pdb --boxcoords -39.355 0.2629999999999999 -4.793 -25.354999999999997 14.263 9.207`


