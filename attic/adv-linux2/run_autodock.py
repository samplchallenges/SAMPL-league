import os
from openmoltools import openeye
from openeye import oechem

PYTHON_PATH = "/opt/app/dependencies/mgl/bin/python"

def SMILES_to_ism(smiles_str: str) -> "filename":
	fname = "/data/smiles.ism"
	smiles_files = open(fname, "w")
	smiles_files.write(smiles_str)
	return fname


def ism_to_pdb(smiles_file: str) -> "filename":
	chgdpdb_file = "/data/chgd_ligand.pdb"
	ifs = oechem.oemolistream()
	if not ifs.open(smiles_file):
		oechem.OEThrow.Fatal("Unable to open %s for reading" % smiles_file)

	for ligand in ifs.GetOEMols():
		ligand = openeye.generate_conformers(ligand)
		ligand = charge.sanitize_OEMol(ligand)
		chgd_ligand = charge.assign_ELF10_charges(ligand)
		
		ofs = oechem.oemolostream()
		if not ofs.open(chgdpdb_file):
			oechem.OEThrow.Fatal("Unable to open %s for reading" % chgdpdb_file)
		ofs.SetFormat(oechem.OEFormat_PDB)
		oechem.OEWriteMolecule(ofs, chgd_ligand)
		ofs.close()
	return chgdpdb_file

def autodock():
	smiles_str = "CCCc1ccccc1"
	ism_file = SMILES_to_ism(smiles_str)
	chgdpdb_file = ism_to_pdb(ism_file)

	
	print("running vina help")
	os.system("/opt/app/dependencies/adv/bin/vina --help")
	print("running prepare_receptor4.py --help")
	os.system(f"{PYTHON_PATH} /opt/app/dependencies/mgl/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l {chgdpdb_file} -o /data/preppedligand.pdbqt")

	os.system(f"{PYTHON_PATH} /opt/app/dependencies/mgl/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r '4w51-cryo.pdb' -o '/data/4w51-cryo.pdbqt'")
	#os.system(f"{python_path} /opt/app/dependencies/mgl/mgltoolspckgs/autodocktools/utilities24/prepare_ligand4.py -l '/data/ligand.pdb'")
	os.system("/opt/app/dependencies/adv/bin/vina --config 'config.txt' --ligand '/data/preppedligand.pdbqt' --out /data/ligand-docked.pdbqt")	


if __name__ == "__main__":
	autodock()
