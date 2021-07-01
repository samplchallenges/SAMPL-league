import mdtraj as md 
from openeye import oechem, oedocking

def SMILES_to_ism(smiles: str) -> "filename":
	fname = "smiles.ism"
	smiles_files = open(fname, "w")
	smiles_files.write(smiles)
	return fname


def PDB_to_oeb(receptor_pdb: str, boxcoords: (float,float,float,float,float,float)) -> "filename":
	receptor_oeb = "/tmp/receptor.oeb"

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
	if boxcoords == None:
		xmin, ymin, zmin, xmax, ymax, zmax = get_xyz_box(receptor_pdb)
	else:
		xmin, ymin, zmin, xmax, ymax, zmax = boxcoords
	print(boxcoords)
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
