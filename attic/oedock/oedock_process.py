# imports for python packages
from openeye import oechem, oedocking
from openmoltools import openeye
import click


# imports for python modules
import charge
import docking
import convertPDB



def open_oemolostream(outstream: str) -> oechem.oemolostream:
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
	if not ofs.open(outstream):
		oechem.OEThrow.Fatal(f"Unable to open {outstream} for writing")
	return ofs



@click.command()
@click.option(
	"-s",
	"--smiles",
	default="CCC",
	help="SMILES string"
)
@click.option(
	"-r",
	"--receptor",
	help="receptor PDB"
)
@click.option(
        "-b",
        "--boxcoords",
        type=click.Tuple([float,float,float,float,float,float]),
        help="The minimum and maximum values of the coordinates of the box representing the binding site.\nEnter in the order: xmin ymin zmin xmax ymax zmax"
)
def dock(smiles: str, receptor: str, boxcoords) -> None:
	receptor_file_path = receptor

	fname = convertPDB.SMILES_to_ism(smiles)

	ifs = oechem.oemolistream()
	if not ifs.open(fname):
		oechem.OEThrow.Fatal("Unable to open %s for reading" % fname)

	for ligand in ifs.GetOEMols():
		ligand = openeye.generate_conformers(ligand)
		ligand = charge.sanitize_OEMol(ligand)
		chgd_ligand = charge.assign_ELF10_charges(ligand)
		print(f"in for loop: {chgd_ligand}")

		receptor_file_path = convertPDB.PDB_to_oeb(receptor, boxcoords)

		dock, sdtag = docking.initialize_docking(receptor_file_path, "Chemgauss4")
		print(f"sdtag: {sdtag}	type: {type(sdtag)}")
		docked_mol, score = docking.dock_molecule(dock, sdtag, ligand)
		print(docked_mol)
		print(f"score: {score}	type: {type(score)}")

		docking.sort_mol_by_score(docked_mol, receptor_file_path)
		ostream = open_oemolostream("dockedmol.oeb.gz")
		oechem.OEWriteMolecule(ostream, ligand)


if __name__ == "__main__":
	"""
	ifs = oechem.oemolistream()
	ifs.SetFormat(oechem.OEFormat_OEB)

	if not ifs.open("sEH-receptor.oeb"):
		oechem.OEThrow.Fatal(f"Unable to open {receptor_pdb} for reading")


	for mol in ifs.GetOEGraphMols():
		print(mol.GetData())
	"""

	dock()
