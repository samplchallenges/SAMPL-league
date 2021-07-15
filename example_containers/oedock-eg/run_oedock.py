# imports for python library
import os

# imports for python packages
from openeye import oechem, oedocking
import openeye.oequacpac as oequacpac
import openeye.oeomega as oeomega
from openeye import oespruce
from openmoltools import openeye
import click

from ligand import Ligand
from receptor import Receptor
from docking import Docking


LIGAND_KEY = "docked_ligand"
RECEPTOR_KEY = "receptor"
SCORE_KEY = "score"


def format_inputs(sz_x,sz_y,sz_z, c_x,c_y,c_z, b_xmin,b_ymin,b_zmin,b_xmax,b_ymax,b_zmax, hint_ligand):
	# Use boxsize and center
	if b_xmin and b_ymin and b_zmin and b_xmax and b_ymax and b_zmax \
			and not sz_x and not sz_y and not sz_z and not c_x and not c_y and not c_z and not hint_ligand:
		boxcoords = (b_xmin, b_ymin, b_zmin, b_xmax, b_ymax, b_zmax)
		hint_ligand = None
		
	elif sz_x and sz_y and sz_z and c_x and c_y and c_z \
			and not b_xmin and not b_ymin and not b_zmin and not b_xmax and not b_ymax and not b_zmax and not hint_ligand:
		boxsize = (sz_x, sz_y, sz_z)
		center = (c_x, c_y, c_z)
		boxcoords = Receptor.get_boxcoords(boxsize, center)
		hint_ligand = None
	elif hint_ligand and not b_xmin and not b_ymin and not b_zmin and not b_xmax and not b_ymax and not b_zmax \
			and not sz_x and not sz_y and not sz_z and not c_x and not c_y and not c_z:
		boxcoords = None
		pass
	else:
		raise Exception("ERROR: boxsize (--sz_x/y/z) and center (--c_x/y/z) OR boxcoords (--b_x/y/zmin and --b_x/y/zmax) OR --hint-ligand must be specified")

	return hint_ligand, boxcoords


@click.command()
@click.option("-s","--smiles",default="CCC",required=True,help="SMILES string")
@click.option("-r","--receptor",required=True,help="receptor PDB")

@click.option("--hint_ligand",help="PDB of ligand docked into receptor to improve docking. Must be used without --boxcoords --boxsize and --center.")

@click.option("--c_x",type=float,help="box center x coordinate; must be used with --c_y/z and --sz_x/y/z")
@click.option("--c_y",type=float,help="box center y coordinate; must be used with --c_x/z and --sz_x/y/z")
@click.option("--c_z",type=float,help="box center z coordinate; must be used with --c_x/y and --sz_x/y/z")

@click.option("--sz_x",type=int,help="box size in the x direction; must be used with --c_x/y/z and --sz_y/z")
@click.option("--sz_y",type=int,help="box size in the y direction; must be used with --c_x/y/z and --sz_x/z")
@click.option("--sz_z",type=int,help="box size in the z direction; must be used with --c_x/y/z and --sz_x/y")

@click.option("--b_xmin",type=float,help="box coordinate x min must be used with --b_ymin/zmin/xmax/ymax/zmax")
@click.option("--b_ymin",type=float,help="box coordinate y min must be used with --b_xmin/zmin/xmax/ymax/zmax")
@click.option("--b_zmin",type=float,help="box coordinate z min must be used with --b_xmin/ymin/xmax/ymax/zmax")
@click.option("--b_xmax",type=float,help="box coordinate x max must be used with --b_xmin/ymin/zmin/ymax/zmax")
@click.option("--b_ymax",type=float,help="box coordinate y max must be used with --b_xmin/ymin/zmin/xmax/zmax")
@click.option("--b_zmax",type=float,help="box coordinate z max must be used with --b_xmin/ymin/zmin/xmax/ymax")

@click.option("--output-dir", help="Output Directory", type=click.Path(exists=True))

def oedock(smiles, receptor, sz_x,sz_y,sz_z, c_x,c_y,c_z, b_xmin,b_ymin,b_zmin,b_xmax,b_ymax,b_zmax, hint_ligand, output_dir) -> None:


	hint_ligand,boxcoords = format_inputs(sz_x,sz_y,sz_z, c_x,c_y,c_z, b_xmin,b_ymin,b_zmin,b_xmax,b_ymax,b_zmax, hint_ligand)

	# make a ligand from smiles string and write to oeb
	ligand_oeb_file = f"{output_dir}/ligand_prepped.oeb"
	lig = Ligand(smiles)
	lig.generate_conformers()
	lig.charge()
	lig.write(ligand_oeb_file, oechem.OEFormat_OEB)


	# make receptor from pdb and write to oeb
	receptor_oeb_file = f"{output_dir}/receptor.oeb"
	rec = Receptor(receptor)
	rec.set_make_receptor_param(hint_ligand, boxcoords)
	rec.make_receptor()
	rec.write(receptor_oeb_file)



	dock_method = oedocking.OEDockMethod_Chemgauss4
	dock_resolution = oedocking.OESearchResolution_Default

	dock = Docking(dock_method, dock_resolution)
	dock.read_receptor_file(receptor_oeb_file)
	dock.initialize()

	num_poses = 1
	dock.dock(ligand_oeb_file,num_poses)


	docked_lig_pdb_file = f"{output_dir}/docked_ligand.pdb"
	score = dock.write_docked_ligands(docked_lig_pdb_file, oechem.OEFormat_PDB)

	docked_rec_pdb_file = f"{output_dir}/docked_receptor.pdb"
	dock.write_receptor(docked_rec_pdb_file, oechem.OEFormat_PDB)

	print(f"{LIGAND_KEY} {docked_lig_pdb_file}")
	print(f"{RECEPTOR_KEY} {docked_rec_pdb_file}")
	print(f"{SCORE_KEY} {score}")

if __name__ == "__main__":
	oedock()
