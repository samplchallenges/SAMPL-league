# Docking Container Requirements

## Input Requirements
* `--receptor`: receptor `.pdb` file to dock the ligand into
  * Example: `--receptor data/receptor.pdb`
* `--smiles`: quoted SMILES string representing the ligand to dock (i.e. "CCC")
  * Example: `--smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O"`
* `--hint`: `.pdb` file with a hint ligand to denote the docking region
  * Example: `--hint data/hint.pdb`
* `--hint_radius`: numeric Angstrom distance from the hint ligand (see above) to consider as the docking region
  * Example: `--hint_radius 6.0`
* `--hint_molinfo`: resname of the hint ligand used in the hint `.pdb` file
  * Example: `--hint_molinfo "E51"`
* `--output-dir`: directory to save final docking files (docked ligand and receptor files)
  * You will not need to handle determining the output directory input as the `ever_given` wrapper handles this for you. Please ensure that the required output files are saved to the output-dir directory

## Output Requirements

**File Outputs**: Output the following files into the `output-dir`
* **docked ligand file**: a `.mol2`, `.pdb` or `.sdf` of the docked ligand
  * path_to_docked_ligand_file = `{output-dir}/{docked_ligand_file}`
* **receptor file**: a `.pdb` file of the receptor used or modified by your docking program, this is important for rmsd scoring purposes in case your complex changes frame of reference
  * path_to_receptor_file = `{output-dir}/{receptor_file}`


**Printed Outputs**: Print the following to to `stdout`
``` docked_ligand {path_to_docked_ligand_file}
receptor {path_to_receptor_file}
```
* Your container should output above in the format `key value` where the keys are `docked_ligand` and `receptor` and the values are file paths. The key and value should be separated by a single space
* These are the only two outputs that should be printed to `stdout`. Please print any extraneous error messages to `stderr` so output parsing is not compromised
* If you are purposely avoiding outputtting a prediction for a compound, please replace `{path_to_docked_ligand_file}` and `{path_to_receptor_file}` with `no_prediction` (see example below)
   ```
   docked_ligand no_prediction
   receptor no_prediction
   ```

## Example Python Main Function Definition
```
import click

@click.command()
@click.option("--receptor", required=True, type=click.Path(exists=True), help="path of receptor PDB to dock the ligand into")
@click.option("--smiles", required=False, type=click.Path(exists=True), help="file with SMILES strings of ligands to be docked")
@click.option("--smiles_arg", required=False, help="SMILES string of a ligand within quotes to avoid command line parsing errors (i.e. \"CCC\")")

@click.option("--hint",required=True,type=click.Path(exists=True),help="path of hint ligand complex for docking region hint")
@click.option("--hint_molinfo",required=True,help="residue name of the ligand in the hint complex")
@click.option("--hint_radius",required=True,type=float,help="box size of the box to dock into")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

def docking_main(receptor, smiles, smiles_argument, hint, hint_molinfo, hint_radius, output_dir):
        ''' docks the given smiles string into the receptor within the area specified by hint and hint-radius
            INPUTS:    receptor:        file    receptor PDB path to dock ligand into
                       smiles:          file    file of SMILES string of ligands to be docked (use either smiles or smiles_arg, not both)
                       smiles_arg:      str     SMILES string of ligand to be docked, use quotes use quotes (use either smiles or smiles_arg, not both)
                       hint:            file    hint PDB contains a receptor ligand complex to show binding site region
                       hint_molinfo:    str     resname of the ligand used in the hint PDB
                       hint_radius:     float   radius around the hint ligand to consider in docking
                       output_dir:      str     output director for receptor and docked_ligand
            OUTPUTS:   prints           docked_ligand {path_to_docked_ligand_file}
                       prints           receptor no_prediction
                       writes file(s)   docked ligand file as a .pdb .mol2 or .sdf
                       writes file(s)   receptor prepped and used by program in docking as .pdb
        '''
```