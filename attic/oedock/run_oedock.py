# imports for python library
import os

# imports for python packages
from openeye import oechem, oedocking
from openmoltools import openeye
import click

# imports for python modules
import charge
import docking
import convertPDB

def get_out_dir(bind_dir: str):
        ''' creates a unique outfile directory to store the results in
        '''
        ID = 0
        while os.path.exists(f"{bind_dir}/out-{ID}"):
                ID += 1
        return f"{bind_dir}/out-{ID}"

def open_oemolostream(outstream: str, oformat) -> oechem.oemolostream:
	''' opens the oechem.oemolostream() with the path outstream
		
	    Parameters:
	    ----------
	    outstream: str
	               String representing the path of the outstream

	    Returns:
	    --------
	    ofs: oechem.oemolostream 
		 outstream to write the molecules to
	'''
	ofs = oechem.oemolostream()
	ofs.SetFormat(oechem.OEFormat_PDB)
	if not ofs.open(outstream):
		oechem.OEThrow.Fatal(f"Unable to open {outstream} for writing")
	return ofs



def get_boxcoords(boxsize: (int,int,int), center: (float,float,float)):
	center_x, center_y, center_z = center
	size_x, size_y, size_z = boxsize

	xmin = center_x - (0.5 * size_x)
	ymin = center_y - (0.5 * size_y)
	zmin = center_z - (0.5 * size_z)
	
	xmax = center_x + (0.5 * size_x)
	ymax = center_y + (0.5 * size_y)
	zmax = center_z + (0.5 * size_z)

	return (xmin, ymin, zmin, xmax, ymax, zmax)	




@click.command()
@click.option(
	"-s",
	"--smiles",
	default="CCC",
	required=True,
	help="SMILES string"
)
@click.option(
	"-r",
	"--receptor",
	required=True,
	help="receptor PDB"
)
@click.option(
        "--boxcoords",
        type=click.Tuple([float,float,float,float,float,float]),
        help="The minimum and maximum values of the coordinates of the box representing the binding site. Enter in the order: xmin ymin zmin xmax ymax zmax"
)
@click.option(
	"--boxsize",
	type=click.Tuple([int,int,int]),
	help="Size of the box to dock into. Enter in the order of: size_x size_y size_z"
)
@click.option(
	"--center",
	type=click.Tuple([float,float,float]),
	help="Center of the box to dock into. Enter in the order of: center_x center_y center_z"
)
@click.option(
	"--bind_out",
	required=True,
	help="directory in the container the outputs are bound to"
)
@click.option(
	"--bind_in",
	required=True,
	help="directory in the container the inputs are bound to"
)
def oedock(smiles: str, receptor: str, boxcoords, boxsize, center, bind_out, bind_in) -> None:

	# Use boxsize and center
	if boxsize != None and center != None and boxcoords == None:
		boxcoords = get_boxcoords(boxsize, center)

	# Use boxcoords
	elif boxsize == None and center == None and boxcoords != None:
       		pass
	
	# Improper inputs
	else: 
		print("ERROR: either --boxcoords alone or the pair of --boxsize and --center together must be provided")
		return

	# Create output directory       
	out_dir = get_out_dir(bind_out)
	os.mkdir(out_dir)

	# Load the ligand smiles as OEMol	
	ligand = oechem.OEMol()
	oechem.OEParseSmiles(ligand, smiles)

	# Assign partial charges and 3D coordinates
	ligand = openeye.generate_conformers(ligand)
	ligand = charge.sanitize_OEMol(ligand)
	chgd_ligand = charge.assign_ELF10_charges(ligand)

	# Turn receptor PDB into oeb
	receptor_file_path = convertPDB.PDB_to_oeb(f"{bind_in}/{receptor}", boxcoords)

	# Dock and score the ligand
	dock, sdtag, receptor = docking.initialize_docking(receptor_file_path, "Chemgauss4")
	docked_mol, score = docking.dock_molecule(dock, sdtag, chgd_ligand)
	print(f"score: {score}	type: {type(score)}")

	# Sort the docked poses by score
	docking.sort_mol_by_score(docked_mol, receptor_file_path)
	
	# Save the docked poses as a PDB
	ofs = open_oemolostream(f"{out_dir}/docked_mol.pdb", oechem.OEFormat_PDB)
	oechem.OEWriteMolecule(ofs, ligand)
	ofs.close()

	# Save the receptor used as a PDB
	ofs = open_oemolostream(f"{out_dir}/docked_receptor.pdb", oechem.OEFormat_PDB)
	oechem.OEWriteMolecule(ofs, receptor)
	ofs.close()

if __name__ == "__main__":
	oedock()
