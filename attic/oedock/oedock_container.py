from openeye import oechem, oedocking
import OpenMMCubes.utils as utils





def initialize_docking(receptor_file_path: str, dock_method_str: str):
	receptor = oechem.OEGraphMol()
	if not oedocking.OEReadReceptorFile( receptor, str( receptor_file ) ):
        	# raise an exception if the receptor file cannot be read
        	raise Exception("Unable to read receptor from {0}".format( receptor_file ))
	dock_method = _set_dock_method( dock_method_str, receptor )
	dock, sdtag = _set_sdtag_and_dock( dock_method )

	if not dock.Initialize(receptor):
        	# raise an exception if the receptor cannot be initialized
		raise Exception("Unable to initialize Docking with {0}".format(self.args.receptor))

	return dock, sdtag


def _set_dock_method( dock_method_str: str, receptor: "OEGraphMol" ) -> "OEDockMethod":
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
    
    if dock_method_str == "Chemgauss4":
        dock_method = oedocking.OEDockMethod_Chemgauss4
        
    elif dock_method_str == "Shapegauss":
        dock_method = oedocking.OEDockMethod_Shapegauss

    elif dock_method_str == "PLP":
        dock_method = oedocking.OEDockMethod_PLP
        
    elif dock_method_str == "Hybrid":
        dock_method = oedocking.OEDockMethod_Hybrid
        if not oedocking.OEReceptorHasBoundLigand(receptor):
            # OEDockMethod_Hybrid requires the receptor to have a ligand bound
            # this option will default to OEDockMethod_Chemgauss4 if no
            # ligand is bound
            oechem.OEThrow.Warning("No bound ligand, switching OEDockMethod to ChemGauss4.")
            dock_method = oedocking.OEDockMethod_Chemgauss4
            
    return dock_method

def _set_sdtag_and_dock( dock_method: "OEDockMethod" ) -> tuple:
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







if __name__ == "__main__":
	dock()
