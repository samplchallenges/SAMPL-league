# python built in
import os
import random

# python packages
from openmoltools import openeye
from openeye import oechem
import click

# my python modules
import charge

PYTHON_PATH = "/opt/app/dependencies/mgl/bin/python"
UTILITIES_PATH = "/opt/app/dependencies/mgl/MGLToolsPckgs/AutoDockTools/Utilities24/"
BIND_DIR = "/data"

ID = random.randint(1000000, 9999999)
INT_DIR = f"{BIND_DIR}/int-{ID}"
OUT_DIR = f"{BIND_DIR}/out-{ID}"


def SMILES_to_ism(smiles_str: str, smiles_path:str):
	smiles_files = open(smiles_path, "w")
	smiles_files.write(smiles_str)


def make_config_file(receptor, smiles, boxsize, center, exhaustiveness, config_path):

	configfile = open(config_path, "w")
	configfile.write(f"receptor = {receptor}\n")

	center_x, center_y, center_z = center
	configfile.write(f"center_x = {center_x}\n")
	configfile.write(f"center_y = {center_y}\n")
	configfile.write(f"center_z = {center_z}\n")

	size_x, size_y, size_z = boxsize
	configfile.write(f"size_x = {size_x}\n")
	configfile.write(f"size_y = {size_y}\n")
	configfile.write(f"size_z = {size_z}\n")

	configfile.write(f"exhaustiveness = {exhaustiveness}")
	configfile.close()




def ism_to_pdb(smiles_path: str, charge_path: str):

	ifs = oechem.oemolistream()
	if not ifs.open(smiles_path):
		oechem.OEThrow.Fatal("Unable to open %s for reading" % smiles_file)

	for ligand in ifs.GetOEMols():
		ligand = openeye.generate_conformers(ligand)
		ligand = charge.sanitize_OEMol(ligand)
		chgd_ligand = charge.assign_ELF10_charges(ligand)
		
		ofs = oechem.oemolostream()
		if not ofs.open(charge_path):
			oechem.OEThrow.Fatal("Unable to open %s for reading" % charge_path)
		ofs.SetFormat(oechem.OEFormat_PDB)
		oechem.OEWriteMolecule(ofs, chgd_ligand)
		ofs.close()

@click.command()
@click.option(
	"-r",
	"--receptor",
	help="path of receptor PDB to dock the ligand into"
)
@click.option(
	"-s",
	"--smiles",
	help="SMILES str of ligand to be docked. quote and add white space at the end \"CCC \""
)
@click.option(
	"-b",
	"--boxsize",
	type=click.Tuple([int,int,int]),
	help="The size of the box to dock into"
)
@click.option(
	"-c",
	"--center",
	type=click.Tuple([float,float,float]),
	help="The center of the box to dock into"
)
@click.option(
	"-e",
	"--exhaustiveness",
	help=""
)

def autodock(receptor, smiles, boxsize, center, exhaustiveness):
	
	os.mkdir(INT_DIR)
	os.mkdir(OUT_DIR)

	ligsmi_path = f"{INT_DIR}/smiles.ism"
	ligchg_path = f"{INT_DIR}/lig-chg.pdbqt"
	ligprep_path = f"{INT_DIR}/lig-prep.pdbqt"

	receptor_path = receptor
	receptorprep_path = f"{INT_DIR}/rec-prep.pdbqt"

	config_path = f"{INT_DIR}/config.txt"

	ligdock_path = f"{OUT_DIR}/lig_dock.pdbqt"


	make_config_file(receptorprep_path, smiles, boxsize, center, exhaustiveness, config_path)

	# SMILES string -> 3D ligand with partial charges assigned -> prepped by autodock
	print("putting smile string in a file")
	SMILES_to_ism(smiles, ligsmi_path)
	print("turning smile into 3d molecule")
	ism_to_pdb(ligsmi_path, ligchg_path)
	print("preparing ligand")
	os.system(f"{PYTHON_PATH} {UTILITIES_PATH}prepare_ligand4.py -l {ligchg_path} -o {ligprep_path}")

	print("preparing receptor")
	# receptor preparation
	os.system(f"{PYTHON_PATH} {UTILITIES_PATH}prepare_receptor4.py -r {receptor_path} -o {receptorprep_path}")


	print("running vina")
	# docking
	print(f"/opt/app/dependencies/adv/bin/vina --config {config_path} --ligand {ligprep_path} --out {ligdock_path}")
	os.system(f"/opt/app/dependencies/adv/bin/vina --config {config_path} --ligand {ligprep_path} --out {ligdock_path}") #\
													#--center_x {center_x} --center_y {center_y} --center_z {center_z} \
													#--size_x {size_x} --size_y {size_y} --size_z {size_z}")


if __name__ == "__main__":
	autodock()
