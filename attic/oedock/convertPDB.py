import mdtraj as md 
from openeye import oechem, oedocking

def SMILES_to_ism(smiles: str) -> "filename":
	fname = "smiles.ism"
	smiles_files = open(fname, "w")
	smiles_files.write(smiles)
	return fname



def get_ifs(filetype: str, filename: str):
	ifs = oechem.oemolistream()
	if filetype == "PDB":
		ifs.SetFormat(oechem.OEFormat_PDB)
	if not ifs.open(filename):
                oechem.OEThrow.Fatal(f"Unable to open {receptor_pdb} for reading")
	return ifs

def PDB_to_mol(pdb: str):
	ifs = get_ifs("PDB", pdb)
	oemol = oechem.OEGraphMol()
	oechem.OEReadMolecule(ifs, oemol)
	ifs.close()
	return oemol



def get_oebox(boxcoords: (float,float,float,float,float,float)):
        if boxcoords == None:
                xmin, ymin, zmin, xmax, ymax, zmax = get_xyz_box(receptor_pdb)
        else:
                xmin, ymin, zmin, xmax, ymax, zmax = boxcoords
        
        return oedocking.OEBox(xmin, ymin, zmin, xmax, ymax, zmax)




def PDB_to_oeb(receptor_pdb: str, ligand_pdb: str, boxcoords: (float,float,float,float,float,float)) -> "filename":
	receptor_oeb_path = "/data/out/receptor.oeb"
	
	protein_oemol = PDB_to_mol(receptor_pdb)
	receptor = oechem.OEGraphMol()

	if ligand_pdb != None:
		ligand_oemol = PDB_to_mol(ligand_pdb)
		oedocking.OEMakeReceptor(receptor, protein_oemol, ligand_oemol)
	else:
		box = get_oebox(boxcoords)
		oedocking.OEMakeReceptor(receptor, protein_oemol, box)
	
	receptor.SetData('MakeReceptor::SiteDetectTag', True)

	ofs = oechem.oemolostream()
	ofs.SetFormat(oechem.OEFormat_OEB)
	if not ofs.open(receptor_oeb_path):
		oechem.OEThrow.Fatal(f"Unable to open {receptor_oeb} for writing")
	oechem.OEWriteMolecule(ofs, receptor)
	ofs.close()
	return receptor_oeb_path

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
