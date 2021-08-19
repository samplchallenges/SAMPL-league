## Input Requirements
* `--receptor`: a receptor `.pdb` file to dock the ligand into
 * Example: `--receptor data/receptor.pdb`
* `--smiles`: a quoted SMILES str representing the ligand to dock (i.e. "CCC")
 * Example: `--smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O"`
* `--hint`: a `.pdb` file with a hint ligand to denote the docking region
 * Example: `--hint data/hint.pdb`
* `--hint_radius`: a numeric Angstrom distance from the hint ligand (see above) to consider as the docking region
 * Example: `--hint_radius 6.0`
* `--hint_molinfo`: the resname of the hint ligand used in the hint `.pdb` file
 * Example: `--hint_molinfo "E51"`
* `--output-dir`: directory to save final docking files (docked ligand and receptor files)
 * You will not need to handle determining the output directory input as the `ever_given` wrapper handles this for you. Please ensure that the required output files are saved to the output-dir directory

## Output Requirements

**File Outputs**: Output the following files into the `output-dir`
* **docked ligand file**: a `.mol2`, `.pdb` or `.sdf` of the docked ligand
  * path_to_docked_ligand_file = `{output-dir}/{docked_ligand_file}`
* **receptor file**: a `.pdb` file of the receptor used or modified by your docking program
  * path_to_receptor_file = `{output-dir}/{receptor_file}`


**Printed Outputs**: Print the following to to `stdout`
```
docked_ligand {path_to_docked_ligand_file}
receptor {path_to_receptor_file}
```
* Your container should output above in the format `key value` where the keys are `docked_ligand`/`receptor` and the values are file paths
* These are the only two outputs that should be printed to stdout. Please print any extraneous error messages to stderr so output parsing is not compromised


