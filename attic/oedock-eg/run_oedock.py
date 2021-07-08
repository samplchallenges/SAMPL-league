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


LIGAND_KEY = "docked_ligand"
RECEPTOR_KEY = "receptor"
SCORE_KEY = "score"


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
	"--hint_ligand",
	help="PDB of ligand docked into receptor to improve docking. Must be used without --boxcoords --boxsize and --center."
)
@click.option("--c_x",type=float)
@click.option("--c_y",type=float)
@click.option("--c_z",type=float)

@click.option("--sz_x",type=int)
@click.option("--sz_y",type=int)
@click.option("--sz_z",type=int)

@click.option("--b_xmin",type=float)
@click.option("--b_ymin",type=float)
@click.option("--b_zmin",type=float)
@click.option("--b_xmax",type=float)
@click.option("--b_ymax",type=float)
@click.option("--b_zmax",type=float)

@click.option("--output-dir", help="Output Directory", type=click.Path(exists=True))

def oedock(smiles, receptor, sz_x,sz_y,sz_z, c_x,c_y,c_z, b_xmin,b_ymin,b_zmin,b_xmax,b_ymax,b_zmax, hint_ligand, output_dir) -> None:

	# Use boxsize and center
	if b_xmin and b_ymin and b_zmin and b_xmax and b_ymax and b_zmax \
			and not sz_x and not sz_y and not sz_z and not c_x and not c_y and not c_z and not hint_ligand:
		boxcoords = (b_xmin, b_ymin, b_zmin, b_xmax, b_ymax, b_zmax)
		hint_ligand = None
		
	elif sz_x and sz_y and sz_z and c_x and c_y and c_z \
			and not b_xmin and not b_ymin and not b_zmin and not b_xmax and not b_ymax and not b_zmax and not hint_ligand:
		boxsize = (sz_x, sz_y, sz_z)
		center = (c_x, c_y, c_z)
		boxcoords = get_boxcoords(boxsize, center)
		hint_ligand = None
	elif hint_ligand and not b_xmin and not b_ymin and not b_zmin and not b_xmax and not b_ymax and not b_zmax \
			and not sz_x and not sz_y and not sz_z and not c_x and not c_y and not c_z:
		boxcoords = None
		pass

	else:
		print("ERROR: boxsize (--sz_x/y/z) and center (--c_x/y/z) OR boxcoords (--b_x/y/zmin and --b_x/y/zmax) OR --hint-ligand must be specified")
		return

	# Load the ligand smiles as OEMol	
	ligand = oechem.OEMol()
	oechem.OEParseSmiles(ligand, smiles)

	# Assign partial charges and 3D coordinates
	ligand = openeye.generate_conformers(ligand)
	ligand = charge.sanitize_OEMol(ligand)
	chgd_ligand = charge.assign_ELF10_charges(ligand)

	# Turn receptor PDB into oeb
	receptor_file_path = convertPDB.PDB_to_oeb(receptor, hint_ligand, boxcoords)
	

	# Dock and score the ligand
	dock, sdtag, receptor = docking.initialize_docking(receptor_file_path, "Chemgauss4")
	docked_mol, score = docking.dock_molecule(dock, sdtag, chgd_ligand)

	# Sort the docked poses by score
	docking.sort_mol_by_score(docked_mol, receptor_file_path)
	
	# Save the docked poses as a PDB
	ofs = open_oemolostream(f"{output_dir}/docked_mol.pdb", oechem.OEFormat_PDB)
	oechem.OEWriteMolecule(ofs, ligand)
	ofs.close()

	# Save the receptor used as a PDB
	ofs = open_oemolostream(f"{output_dir}/docked_receptor.pdb", oechem.OEFormat_PDB)
	oechem.OEWriteMolecule(ofs, receptor)
	ofs.close()

	print(f"{LIGAND_KEY} {output_dir}/docked_mol.pdb")
	print(f"{RECEPTOR_KEY} {output_dir}/docked_receptor.pdb")
	print(f"{SCORE_KEY} {score}", end="")

if __name__ == "__main__":
	oedock()
