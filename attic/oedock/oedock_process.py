import charge
import docking
from openeye import oechem, oedocking
from openmoltools import openeye
import click
import mdtraj as md 


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


def SMILES_to_ism(smiles: str) -> "filename":
	fname = "smiles.ism"
	smiles_files = open(fname, "w")
	smiles_files.write(smiles)
	return fname


def PDB_to_oeb(receptor_pdb: str) -> "filename":
	receptor_oeb = "receptor.oeb"

	ifs = oechem.oemolistream()
	ofs = oechem.oemolostream()

	ifs.SetFormat(oechem.OEFormat_PDB)
	ofs.SetFormat(oechem.OEFormat_OEB)

	if not ifs.open(receptor_pdb):
		oechem.OEThrow.Fatal(f"Unable to open {receptor_pdb} for reading")
	if not ofs.open(receptor_oeb):
		oechem.OEThrow.Fatal(f"Unable to open {receptor_oeb} for reading")


	# from: https://github.com/choderalab/kinase-resistance-mutants/blob/master/docking/docking.py
	protein_oemol = oechem.OEGraphMol()
	oechem.OEReadMolecule(ifs, protein_oemol)

	xmin, ymin, zmin, xmax, ymax, zmax = get_xyz_box(receptor_pdb)
	box = oedocking.OEBox(xmin, ymin, zmin, xmax, ymax, zmax)

	receptor = oechem.OEGraphMol()

	oedocking.OEMakeReceptor(receptor, protein_oemol, box)

	receptor.SetData('MakeReceptor::SiteDetectTag', True)
	oechem.OEWriteMolecule(ofs, receptor)
	print(receptor.GetData())

	return receptor_oeb

def get_xyz_box(receptor_pdb:str):
	traj=md.load_pdb(receptor_pdb)

	xcoords = traj.xyz[0][:,[0]]
	ycoords = traj.xyz[0][:,[1]]
	zcoords = traj.xyz[0][:,[2]]

	xmin = float(min(xcoords)[0])
	xmax = float(max(xcoords)[0])
	ymin = float(min(ycoords)[0])
	ymax = float(max(ycoords)[0])
	zmin = float(min(zcoords)[0])
	zmax = float(max(zcoords)[0])

	return xmin, ymin, zmin, xmax, ymax, zmax

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
def dock(smiles: str, receptor: str) -> None:
	receptor_file_path = receptor

	fname = SMILES_to_ism(smiles)

	ifs = oechem.oemolistream()
	if not ifs.open(fname):
		oechem.OEThrow.Fatal("Unable to open %s for reading" % fname)

	for ligand in ifs.GetOEMols():
		ligand = openeye.generate_conformers(ligand)
		ligand = charge.sanitize_OEMol(ligand)
		chgd_ligand = charge.assign_ELF10_charges(ligand)
		print(f"in for loop: {chgd_ligand}")

		receptor_file_path = PDB_to_oeb(receptor)

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