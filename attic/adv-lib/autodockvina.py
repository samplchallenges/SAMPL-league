import click
from vina import Vina
import os
from os import path

def file_exists(filename, filetype):
	if not filename:
		raise FileNotFoundError(f"No {filetype} file provided")
	
	if not path.isfile(filename):
		raise FileNotFoundError(f"File \"{filename}\" does not exist")



@click.command()
@click.option(
	"-r", 
	"--receptor", 
	type=str,
	help="receptor pdbqt file" 
)
@click.option(
	"-l", 
	"--ligand", 
	type=str,
	help="ligand pdbqt file"
)
@click.option(
	"-c", 
	"--center", 
	type=click.Tuple([float,float,float]),
	default=(0.,0.,0.),
	help="center of the binding site"
)
@click.option(
	"-b", 
	"--boxsize", 
	type=click.Tuple([int,int,int]),
	default=(30,30,30),
	help="box size of the binding site area" 
)
@click.option(
	"-o",
	"--outfile",
	default="docking-results.pdbqt",
	help="name of docking results pdbqt file"
)
def run_autodock(receptor, ligand, center, boxsize, outfile):
	print("receptor file: ", receptor)
	print("lig file:      ", ligand)
	print("center:        ", center)
	print("boxsize:       ", boxsize)
	print("outfile:       ", outfile)

	try: 
		file_exists(receptor, "receptor")
		file_exists(ligand, "ligand")
	except FileNotFoundError as error:
		print(f"ERROR: {error}")
		return
	
	vna = Vina()
	vna.set_receptor(rigid_pdbqt_filename=receptor)
	vna.set_ligand_from_file(ligand)
	vna.compute_vina_maps(center=list(center), box_size=list(boxsize))
	print(vna.score())
	print(vna.optimize())
	vna.dock(exhaustiveness=1)
	#lig_basename = ligand[:-7]
	vna.write_poses(pdbqt_filename=outfile)


if __name__ == "__main__":
	run_autodock()
