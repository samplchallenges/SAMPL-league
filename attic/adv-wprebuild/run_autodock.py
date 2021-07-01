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


def get_out_dir(bind_dir: str):
	''' creates a unique outfile directory to store the results in
	'''
	ID = 0
	while os.path.exists(f"{bind_dir}/out-{ID}"):
		ID += 1
	return f"{bind_dir}/out-{ID}"


def smiles_to_smi(smiles_str: str, smiles_path:str):
	''' transfer SMILES string to a .smi file
	    used for debugging
	'''
	smiles_files = open(smiles_path, "w")
	smiles_files.write(smiles_str)


def make_config_file(receptor, smiles, boxsize, center, exhaustiveness, config_path):
	''' makes the configuration file which is an input to the vina 
	    command
	'''
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
	''' converts a .smi file to a pdb of a charged ligand with 3D coordinates
	    pairs with smiles_to_smi(smiles, smiles_path) for debugging purposes
	'''
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


def smiles_to_chgpdb(smiles: str, charge_path: str):
	''' convers a smiles string to a pdb of a charged ligand with 3D coordinates
	'''
	ligand = oechem.OEMol()
	oechem.OEParseSmiles(ligand, smiles)
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
	''' gets the last atom index used the receptor pdb to append 
	    ligand to receptor pdb (no longer necessary)
	'''
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
	''' appends the ligand receptor with correct indexes to the end of the receptor pdb
	    (no longer necessary)
	'''
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


def pdbqt_to_pdb(pdbqt_path, pdb_path):
	os.system(f"cut -c-66 {pdbqt_path} > {pdb_path}")


def clean_tmp():
	os.system("rm /tmp/*")

def get_center_and_boxsize(boxcoords):
	xmin, ymin, zmin, xmax, ymax, zmax = boxcoords

	center_x = 0.5 * (xmin + xmax)
	center_y = 0.5 * (ymin + ymax)
	center_z = 0.5 * (zmin + zmax)

	size_x = xmax - xmin
	size_y = ymax - ymin
	size_z = zmax - zmin

	return ((center_x, center_y, center_z), (size_x, size_y, size_z))


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
	"--boxcoords",
	type=click.Tuple([float,float,float,float,float,float])
)
@click.option(
	"-e",
	"--exhaustiveness",
	help=""
)
@click.option(
	"--bind_out",
	help="directory the outputs should be bound to"
)
@click.option(
	"--bind_in",
	help="directory the inputs should be bound to"
)
def autodock(receptor, smiles, boxsize, center, boxcoords, exhaustiveness, bind_out, bind_in):
	''' docks the given smiles string into the receptor within the box specified by
	    boxsize and center
	    exhaustiveness does not work at this moment
	'''

	if boxcoords != None and center != None or boxcoords != None and boxsize != None:
		print("ERROR: either --boxcoords alone or the pair --center and --boxsize together can be used")
		return
	if center != None and boxsize == None or center == None and boxsize != None:
		print("ERROR: --center and --boxsize must be used together")
		return
	if boxcoords != None:
		center, boxsize = get_center_and_boxsize(boxcoords)



	# Create output directory	
	out_dir = get_out_dir(bind_out)
	os.mkdir(out_dir)

	# Set file names for intermediate and output files
	# Output file names should probably become an input option
	ligchg_sdf_path = "/tmp/lig-chg.sdf"
	ligchg_pdbqt_path = "/tmp/lig-chg.pdbqt"
	ligprep_path = "/tmp/lig-prep.pdbqt"

	receptor_path = f"{bind_in}/{receptor}"
	receptorprep_pdbqt_path = "/tmp/rec-prep.pdbqt"
	receptor_pdb_path = f"{out_dir}/rec-dock.pdb"

	config_path = f"{out_dir}/config.txt"

	ligdock_path = f"{out_dir}/lig_dock.pdbqt"
	highscore_pdbqt_path = "/tmp/best_dock.pdbqt"
	highscore_pdb_path = f"{out_dir}/best_dock.pdb"

	complex_pdb_path = f"{out_dir}/docked_complex.pdb"

	# Make the config file
	make_config_file(receptorprep_pdbqt_path, smiles, boxsize, center, exhaustiveness, config_path)

	
	# SMILES string -> 3D ligand with partial charges assigned -> autodock prepped ligand
	print("LOAD: loading smiles string and converting to 3D sdf")
	smiles_to_chgpdb(smiles, ligchg_sdf_path)
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

	# Save the highest scoring pose which should be the first pose in the file
	# and convert this pose to a pdb
	print("SAVE: saving highest scoring pose")
	save_highest_score(ligdock_path, highscore_pdbqt_path)
	pdbqt_to_pdb(highscore_pdbqt_path, highscore_pdb_path)	

	# Convert prepped receptor from pdbqt to pdb
	print("SAVE: saving protein ligand complex")
	pdbqt_to_pdb(receptorprep_pdbqt_path, receptor_pdb_path)
	
	# Remove intermediate files
	print("CLEAN: clearing tmp files")
	clean_tmp()
	
if __name__ == "__main__":
	autodock()
