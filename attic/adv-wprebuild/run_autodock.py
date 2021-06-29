# python built in
import os
import random

# python packages
from openmoltools import openeye
from openeye import oechem
import click
import mdtraj as md

# my python modules
import charge

PYTHON_PATH = "/opt/app/dependencies/mgl/bin/python"
UTILITIES_PATH = "/opt/app/dependencies/mgl/MGLToolsPckgs/AutoDockTools/Utilities24/"
VINA_PATH = "/opt/app/dependencies/adv/bin/vina"
BIND_DIR = "/data"

ID = random.randint(1000000, 9999999)
INT_DIR = f"{BIND_DIR}/int-{ID}"
OUT_DIR = f"{BIND_DIR}/out-{ID}"


def SMILES_to_smi(smiles_str: str, smiles_path:str):
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

	#configfile.write(f"exhaustiveness = {exhaustiveness}")
	configfile.close()




def smi_to_pdb(smiles_path: str, charge_path: str):

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
		ofs.SetFormat(oechem.OEFormat_SDF)
		oechem.OEWriteMolecule(ofs, chgd_ligand)
		ofs.close()


def save_highest_score(dock_path: str, highscore_path: str):
	''' writes the highest scoring docked pose into its own file
	'''

	infile = open(dock_path, 'r')
	outfile = open(highscore_path, 'w')
	
	for line in infile.readlines():
		outfile.write(line)
		if "ENDMDL" in line:
			break
			
	outfile.close()
	infile.close()


def get_last_receptor_ind(receptor_path: str):
	with open(receptor_path) as f:
		for line in f:
			pass
		last_line = line
		splt = last_line.split()
		if splt[0] == "TER":
			return int(splt[1]) - 1
		else:
			return int(split[1])

def append_ligand_receptor(complex_path: str, receptor_path: str, ligand_path: str, lastind: int):
	complexf = open(complex_path, 'w')
	with open(receptor_path) as receptorf:
		for line in receptorf:
			if "TER" in line[0:3]:
				complexf.write("TER\n")
			else:
				complexf.write(line)
	with open(ligand_path) as ligandf:
		ct = 1
		for line in ligandf:
			if "ATOM" in line[0:4]:
				splt = line.split()
				for i in range(len(splt)):
					if i == 1:
						complexf.write(str(lastind + ct) + "\t")
					else:
						complexf.write(splt[i] + "\t")	
				complexf.write("\n")
				ct += 1
		complexf.write("END")
	complexf.close()
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

	ligsmi_path = f"{INT_DIR}/smiles.smi"
	ligchg_sdf_path = f"{INT_DIR}/lig-chg.sdf"
	ligchg_pdbqt_path = f"{INT_DIR}/lig-chg.pdbqt"
	ligprep_path = f"{INT_DIR}/lig-prep.pdbqt"

	receptor_path = receptor
	receptorprep_pdbqt_path = f"{INT_DIR}/rec-prep.pdbqt"
	receptor_pdb_path = f"{OUT_DIR}/rec-dock.pdb"

	config_path = f"{INT_DIR}/config.txt"

	ligdock_path = f"{OUT_DIR}/lig_dock.pdbqt"
	highscore_pdbqt_path = f"{OUT_DIR}/best_dock.pdbqt"
	highscore_pdb_path = f"{OUT_DIR}/best_dock.pdb"

	complex_pdb_path = f"{OUT_DIR}/docked_complex.pdb"

	make_config_file(receptorprep_pdbqt_path, smiles, boxsize, center, exhaustiveness, config_path)

	# SMILES string -> 3D ligand with partial charges assigned -> prepped by autodock
	print("LOAD: putting smile string in a file")
	SMILES_to_smi(smiles, ligsmi_path)
	print("PREP: adding partial charges and 3d coords")
	smi_to_pdb(ligsmi_path, ligchg_sdf_path)
	print("PREP: converting sdf to pdbqt using openbabel")
	os.system(f"obabel {ligchg_sdf_path} -O {ligchg_pdbqt_path}")
	print("PREP: preparing ligand")
	os.system(f"{PYTHON_PATH} {UTILITIES_PATH}prepare_ligand4.py -l {ligchg_pdbqt_path} -o {ligprep_path}")
	

	# Preparing receptor
	print("PREP: preparing receptor")
	os.system(f"{PYTHON_PATH} {UTILITIES_PATH}prepare_receptor4.py -r {receptor_path} -o {receptorprep_pdbqt_path}")


	# Running AutoDock Vina
	print("DOCK: running vina")
	print(f"{VINA_PATH} --config {config_path} --ligand {ligprep_path} --out {ligdock_path}")
	os.system(f"{VINA_PATH} --config {config_path} --ligand {ligprep_path} --out {ligdock_path}")

	
	print("SAVE: saving highest scoring pose")
	# Save the highest scoring pose which should be the first pose in the file
	save_highest_score(ligdock_path, highscore_pdbqt_path)
	os.system(f"cut -c-66 {highscore_pdbqt_path} > {highscore_pdb_path}")	

	print("SAVE: saving protein ligand complex")
	# convert prepped receptor from pdbqt to pdb
	os.system(f"cut -c-66 {receptorprep_pdbqt_path} > {receptor_pdb_path}")

	#os.system(f"cat {receptor_pdb_path} {highscore_pdb_path} | grep -v '^END   ' | grep -v '^END$' >  {complex_pdb_path}")

	#lastind = get_last_receptor_ind(receptor_pdb_path)

	#append_ligand_receptor(complex_pdb_path, receptor_pdb_path, highscore_pdb_path, lastind)


if __name__ == "__main__":
	autodock()
