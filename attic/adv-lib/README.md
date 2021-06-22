To build:
`docker build -t vina .`


To run:
`docker run -it --rm -v $(pwd):/data vina -r /data/<receptor_path> -l /data/<ligand_path>`


