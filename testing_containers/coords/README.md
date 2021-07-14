# How to build


Run `docker build -t calc-coords:0.1 .` to build the container.
Replace if you want to push the container to the hub, include your username in the tag (for example if your user name is mmh42) you would run `docker build -t mmh42/calc-coords:0.1 .`

# How to use

Run `docker run -it --rm calc-coords:0.1 CCC` to calculate the 3D coords of propane (CCC).
You can calculate the molecular mass of other molecules by using their respective smiles string.


Or

docker run --mount 'type=bind,src=/Users/robbason/mobley/SAMPL-league/example_containers/coords/dfile,dst=/mnt/inputfiles,readonly' --mount 'type=bind,src=/Users/robbason/mobley/SAMPL-league/example_containers/output,dst=/mnt/outputfiles' -it robbason/calc-coords:latest --molfile /mnt/inputfiles/ChEBI_16716.mdl --output-dir /mnt/outputfiles

# More information

See the documentation here: https://docs.docker.com/engine/reference/commandline/run/ for more information on the `docker run` command line options and see the documentation here: https://docs.docker.com/engine/reference/commandline/build/ for more information about the `docker build` command.
