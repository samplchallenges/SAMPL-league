# python built in
import os
import random
import sys

# python packages
import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem

SCORE_KEY = "score"
LIGAND_KEY = "docked_ligand"
RECEPTOR_KEY = "receptor"


PYTHON_PATH = "/opt/app/dependencies/mgl/bin/python"
UTILITIES_PATH = "/opt/app/dependencies/mgl/MGLToolsPckgs/AutoDockTools/Utilities24/"
VINA_PATH = "/opt/app/dependencies/adv/bin/vina"


def print_debug(debug: bool, msg:str):
	print(f"{msg}\n" if debug else "", end="")
	sys.stdout.flush()

def get_out_dir(bind_dir: str):
	''' creates a unique outfile directory to store the results in
	'''
	ID = 0
	while os.path.exists(f"{bind_dir}/out-{ID}"):
		ID += 1
	return f"{bind_dir}/out-{ID}"



def make_config_file(receptor, smiles, flex, boxsize, center, exhaustiveness, num_modes, config_path):
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

	if exhaustiveness != None:
		configfile.write(f"exhaustiveness = {exhaustiveness}\n")
	if num_modes != None:
		configfile.write(f"num_modes = {num_modes}")
	
	configfile.close()



def smiles_to_chgpdb(smiles: str, charge_path: str):
	''' convers a smiles string to a pdb of a charged ligand with 3D coordinates
	'''

	mol = Chem.MolFromSmiles(smiles)
	mol2 = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol2)

	ostream = Chem.SDWriter(charge_path)
	ostream.write(mol2)


def save_highest_score(dock_path: str, highscore_path: str):
	''' writes the highest scoring docked pose into its own file
	'''

	infile = open(dock_path, 'r')
	outfile = open(highscore_path, 'w')
	
	for line in infile.readlines():
		if "MODEL" in line:
			continue
		if "ENDMDL" in line:
			break
		outfile.write(line)
			
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

def get_score(score_path):
	with open(score_path) as scoref:
		for line in scoref:
			if "Aff" in line:
				return float(line.split()[1])


def get_utility_cmd(pyfile: str):
	return f"{PYTHON_PATH} {UTILITIES_PATH}{pyfile}"


@click.command()
@click.option("-r","--receptor",help="path of receptor PDB to dock the ligand into")
@click.option("-s","--smiles",help="SMILES str of ligand to be docked. quote and add white space at the end \"CCC \"")
@click.option("--flex",help="Not working yet. flexible sidechains if any pdb")
@click.option("--output-dir", help="Output Directory", type=click.Path(exists=True))

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

@click.option("-n","--num_modes",type=int,help="Number of modes to dock")
@click.option("-e","--exhaustiveness",type=int,help="exhaustiveness of the global search, default=8")

@click.option('--debug', is_flag=True,help="prints debug print statements when --debug flag is used")

def autodock(receptor, smiles, flex, sz_x,sz_y,sz_z, c_x,c_y,c_z, b_xmin,b_ymin,b_zmin,b_xmax,b_ymax,b_zmax, exhaustiveness, num_modes, output_dir, debug):
	''' docks the given smiles string into the receptor within the box specified by
	    boxsize and center
	    exhaustiveness does not work at this moment
	'''

	# Set file names for intermediate and output files
	# Output file names should probably become an input option
	ligchg_sdf_path = "/tmp/lig-chg.sdf"
	ligchg_pdbqt_path = "/tmp/lig-chg.pdbqt"
	ligprep_path = "/tmp/lig-prep.pdbqt"

	receptor_path = f"{receptor}"
	receptorprep_pdbqt_path = "/tmp/rec-prep.pdbqt"
	receptor_pdb_path = f"{output_dir}/rec-dock.pdb" # switch to join()

	config_path = "/tmp/config.txt"

	ligdock_path = "/tmp/lig_dock.pdbqt"
	highscore_pdbqt_path = "/tmp/best_dock.pdbqt"
	highscore_pdb_path = f"{output_dir}/best_dock.pdb" # switch to join()

	score_path = "/tmp/score.txt"

	# Make the config file
	if b_xmin and b_ymin and b_zmin and b_xmax and b_ymax and b_zmax:
		boxcoords = (b_xmin, b_ymin, b_zmin, b_xmax, b_ymax, b_zmax)
		boxsize,center = get_center_and_boxsize(boxcoords)

	elif sz_x and sz_y and sz_z and c_x and c_y and c_z:
		boxsize = (sz_x, sz_y, sz_z)
		center = (c_x, c_y, c_z)
	
	else:
		print("ERROR: Either boxsize (sz_x/y/z) and center (c_x/y/z) OR boxcoords (b_x/y/zmin and b_x/y/zmax) must be specified")

	make_config_file(receptorprep_pdbqt_path, smiles, flex, boxsize, center, exhaustiveness, num_modes, config_path)

	
	# SMILES string -> 3D ligand with partial charges assigned -> autodock prepped ligand
	print("LOAD:   loading smiles string and converting to 3D sdf\n" if debug else "", end="")
	smiles_to_chgpdb(smiles, ligchg_sdf_path)
	print("PREP:   converting sdf to pdbqt using openbabel\n" if debug else "", end="")
	os.system(f"obabel {ligchg_sdf_path} -O {ligchg_pdbqt_path} 2>/dev/null 1>/dev/null")
	print("PREP:   preparing ligand\n" if debug else "", end="")
	cmd = get_utility_cmd("prepare_ligand4.py")
	os.system(f"{cmd} -l {ligchg_pdbqt_path} -o {ligprep_path} > /dev/null")
	

	# Preparing receptor
	print_debug(debug, "PREP:   preparing receptor")
	cmd = get_utility_cmd("prepare_receptor4.py")
	os.system(f"{cmd} -r {receptor_path} -o {receptorprep_pdbqt_path} > /dev/null")


	# Running AutoDock Vina
	print_debug(debug, "DOCK:   running vina")
	#print(f"{VINA_PATH} --config {config_path} --ligand {ligprep_path} --out {ligdock_path}")
	os.system(f"{VINA_PATH} --config {config_path} --ligand {ligprep_path} --out {ligdock_path} > /dev/null")
	
	# Save the highest scoring pose which should be the first pose in the file
	# and convert this pose to a pdb
	print_debug(debug, "SAVE:   saving highest scoring pose")
	save_highest_score(ligdock_path, highscore_pdbqt_path)
	pdbqt_to_pdb(highscore_pdbqt_path, highscore_pdb_path)	

	# Get score with more sigfigs for highest scoring pose
	print_debug(debug, "SCORE:  scoring highest scorint pose")
	os.system(f"{VINA_PATH} --config {config_path} --ligand {highscore_pdbqt_path} --out /tmp/rescore.pdbqt --score_only | grep 'Affinity:' > {score_path}")
	score = get_score(score_path)

	# Convert prepped receptor from pdbqt to pdb
	print_debug(debug, "SAVE:   saving protein ligand complex")
	pdbqt_to_pdb(receptorprep_pdbqt_path, receptor_pdb_path)
	
	# Remove intermediate files
	print_debug(debug, "CLEAN:  clearing tmp files")
	clean_tmp()

	print(f"{SCORE_KEY} {score} kcal/mol")
	print(f"{LIGAND_KEY} {highscore_pdb_path}")
	print(f"{RECEPTOR_KEY} {receptor_pdb_path}", end="")
if __name__ == "__main__":
	autodock()
