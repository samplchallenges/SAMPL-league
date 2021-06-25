from openeye import oechem, oedocking, oequacpac
from openmoltools import openeye

###################################
#### Ligand Charging Functions ####
###################################


def assign_ELF10_charges(mol: oechem.OEMol, max_confs:int=800, strictStereo:bool=True) -> oechem.OEMol:
	''' assigns partial charges to a copy of an uncharged OEMol object
		This function was pulled from: 
			nathanmlim/blues-apps/OrionFloes/LigPrepCubes/ff_utils.py > assignELF10charges fx
	
		Parameters:
		-----------
		mol: oechem.OEMol  
			The molecule that needs to be chaged
		max_confs: int
			The max number of conformers used to calculate the atomic partial charges; default=800
		strict_stereo: bool
			Flag used to check if atoms need to have assigned stereo chemistry or not; default=True

		Returns:
		--------
		mol_copy: oechem.OEMol
			A copy of the original molecule with assigned atomic partial charges
	'''
	mol_copy = mol.CreateCopy()

	if not mol_copy.GetMaxConfIdx() > 200:
		# Generate up to max_confs conformers
 		mol_copy = openeye.generate_conformers(mol_copy, max_confs=max_confs, strictStereo=strictStereo)

	# Assign MMFF Atom types
	if not oechem.OEMMFFAtomTypes(mol_copy):
		raise RuntimeError("MMFF atom type assignment returned errors")

	# ELF10 charges
	status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCELF10Charges())

	if not status:
		raise RuntimeError("OEAssignCharges returned error code %d" % status)

	return mol_copy

def sanitize_OEMol(mol: oechem.OEMol) -> oechem.OEMol:
	''' This function checks if the molecule has coordinates,
		explicit hydrogens and aromaticity. If the molecule
		does not have coordinates a fatal error is raised.
		If the molecule does not have hydrogens or aramatic
		flags are missing then a copy of the molecule is fixed
		Parameters:
		-----------
		mol: oechem.OEMol
			The molecule to be checked
		Return:
		-------
		mol_copy: oechem.OEMol
			A copy of the checked molecule with fixed aromaticity and hydrogens
	'''
	mol_copy = mol.CreateCopy()

	# Check if the molecule has 3D coordinates
	if not oechem.OEGetDimensionFromCoords(mol_copy):
		oechem.OEThrow.Fatal("The molecule coordinates are set to zero")
	# Check if the molecule has hydrogens
	if not oechem.OEHasExplicitHydrogens(mol_copy):
		oechem.OEAddExplicitHydrogens(mol_copy)
	# Check if the molecule has assigned aromaticity
	if not mol_copy.HasPerceived(oechem.OEPerceived_Aromaticity):
		oechem.OEAssignAromaticFlags(mol_copy, oechem.OEAroModelOpenEye)

	# TEMPORARY PATCH FOR SMIRNOFF
	oechem.OETriposAtomNames(mol_copy)
	# Check for any missing atom names, if found reassign all of them
	# if any([atom.GetName() == '' for atom in mol_copy.GetAtoms()]):
	#     oechem.OETriposAtomNames(mol_copy)

	return mol_copy