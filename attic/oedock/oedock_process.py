import sys
from openeye import oechem, oedocking, oequacpac
from openmoltools import openeye


###################################
#### Ligand Charging Functions ####
###################################

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



def smiles_to_OEMol(smiles: str) -> oechem.OEMol:
	''' converts a SMILES string to an OEMol object

	'''
	mol = oechem.OEMol()
	oechem.OESmilesToMol(mol, smiles)
	return mol


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

  
###################################
#### Ligand Docking Functions ####
###################################


def initialize_docking( receptor_file_path: str, dock_method_str: str ) -> tuple:
	''' Sets the receptor, Initializes the receptor for docking and returns the
		dock and sdtag as a tuple i.e. ( dock, sdtag )
		parameters:
			- receptor_file_path: string for the file path to the receptor
			- dock_method_str: string for the dock method to be used
	'''
	receptor = oechem.OEGraphMol()

	if not oedocking.OEReadReceptorFile( receptor, str( receptor_file_path ) ):
		# raise an exception if the receptor file cannot be read
		raise Exception("Unable to read receptor from {0}".format( receptor_file_path ))

	dock_method = set_dock_method( dock_method_str, receptor )

	dock, sdtag = set_sdtag_and_dock( dock_method )

	if not dock.Initialize(receptor):
		# raise an exception if the receptor cannot be initialized
		raise Exception("Unable to initialize Docking with {0}".format(self.args.receptor))

	return dock, sdtag

def dock_molecule( dock: "OEDock", sdtag: str, mcmol, num_poses: int=10 ) -> (oechem.OEMol, float) :
	''' Docks the multiconfomer molecule, with the given number of poses
		Returns a tuple of the docked molecule (dockedMol) and its score
		i.e. ( dockedMol, score )
		parameters:
			- dock: OEDock object
			- sdtag: string representing the name of the docking method
			- numpose: int with the number of poses to generate for each ligand
			- mcmol: multicomformer molecule
	'''
	dockedMol = oechem.OEMol()

	# dock, annotate, and score the molecule
	res = dock.DockMultiConformerMolecule(dockedMol, mcmol, num_poses)
	if res == oedocking.OEDockingReturnCode_Success:
		oedocking.OESetSDScore(dockedMol, dock, sdtag)
		dock.AnnotatePose(dockedMol)
		score = dock.ScoreLigand(dockedMol)
		oechem.OESetSDData(dockedMol, sdtag, "{}".format(score))
		_clean(dockedMol)
		return dockedMol, score
    
	else:
		# raise an exception if the docking is not successful
		raise Exception("Unable to dock ligand {0} to receptor".format( dockedMol ))


def set_dock_method( dock_method_str: str, receptor: oechem.OEGraphMol ) -> "oedock.OEDockMethod":
	''' sets the dock method based on the input parameter "dock_method_str"
		OPTIONS: "Chemgauss4" = OEDockMethod_Chemgauss4
				 "Shapegauss" = OEDockMethod_Shapegauss
				 "PLP" = oedocking.OEDockMethod_PLP
				 "Hybrid" = OEDockMethod_Hybrid
		the dock method may be switched if the requirements for the docking method
		are not met (i.e. Hybrid switches to Chemgauss4 when no ligand is bound)
		returns the dock method
		Parameters
		- dock_method_str: string specifying the dock method to be used
		- receptor: OEGraphMol object
	'''

	dock_method_dict = {
		"Chemgauss4": oedocking.OEDockMethod_Chemgauss4,
		"Shapegauss": oedocking.OEDockMethod_Shapegauss,
		"PLP": oedocking.OEDockMethod_PLP,
		"Hybrid": oedocking.OEDockMethod_Hybrid
	}
        
	if dock_method_str == "Hybrid":
		dock_method = dock_method_dict["Hybrid"]
		if not oedocking.OEReceptorHasBoundLigand(receptor):
			# OEDockMethod_Hybrid requires the receptor to have a ligand bound
			# this option will default to OEDockMethod_Chemgauss4 if no
 			# ligand is bound
			oechem.OEThrow.Warning("No bound ligand, switching OEDockMethod to ChemGauss4.")
			dock_method = dock_method_dict["Chemgauss4"]

	else:
		try:
			dock_method = dock_method_dict[dock_method_str]
		except KeyError:
			oechem.OEThrow.Fatal("Invalid docking method.")
	return dock_method

def set_sdtag_and_dock( dock_method: "oedocking.OEDockMethod" ) -> (oedocking.OEDock, str):
	''' sets the dock method, then sets the sdtag and dock
		Returns a tuple with dock as the first item and sdtag as the second
		i.e. ( dock, sdtag )
		Parameters:
			- dock_method: OEDockMethod object
				- Ex: OEDockMethod_Chemgauss4, OEDockMethod_Shapegauss, 
					  OEDockMethod_PLP, OEDockMethod_Hybrid
	'''
	dock_resolution = oedocking.OESearchResolution_Default
	sdtag = oedocking.OEDockMethodGetName( dock_method )
	dock = oedocking.OEDock( dock_method, dock_resolution)
	return dock, sdtag

def sort_mol_by_score( dockedMol: oechem.OEMol, receptor_file_path: str) -> None:
	''' sorts the poses for each molecule by their score, where higher scores are better
		Parameters:
			- dockedMol: OEMol object
			- receptor_file_path: string file path to a receptor file
	'''
	receptor = oechem.OEGraphMol( )
	oedocking.OEReadReceptorFile( receptor, receptor_file_path )
	score = oedocking.OEScore(oedocking.OEScoreType_Chemgauss4)
	score.Initialize(receptor)
	score.SystematicSolidBodyOptimize(dockedMol)
	score.AnnotatePose(dockedMol)
	sdtag = score.GetName()
	oedocking.OESetSDScore(dockedMol, score, sdtag)
	oechem.OESortConfsBySDTag(dockedMol, sdtag, score.GetHighScoresAreBetter())
	return None


def _clean(mol: "OEMol") -> None:
	''' Parameters:
            - mol: OEMol object
	'''
	mol.DeleteData('CLASH')
	mol.DeleteData('CLASHTYPE')
	mol.GetActive().DeleteData('CLASH')
	mol.GetActive().DeleteData('CLASHTYPE')


def dock():
	receptor_file_path = "sEH-receptor.oeb"

	ifs = oechem.oemolistream()
	if not ifs.open("6NX1.ism"):
		oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[1])

	for ligand in ifs.GetOEMols():
		ligand = openeye.generate_conformers(ligand)
		ligand = sanitize_OEMol(ligand)
		chgd_ligand = assign_ELF10_charges(ligand)
		print(f"in for loop: {chgd_ligand}")
		dock, sdtag = initialize_docking(receptor_file_path, "Chemgauss4")
		print(f"sdtag: {sdtag}	type: {type(sdtag)}")
		docked_mol, score = dock_molecule(dock, sdtag, ligand)
		print(docked_mol)
		print(f"score: {score}	type: {type(score)}")

		sort_mol_by_score(docked_mol, receptor_file_path)
		ostream = open_oemolostream("dockedmol.oeb.gz")
		oechem.OEWriteMolecule(ostream, ligand)
