## Input Requirements
* `--receptor`: a receptor `.pdb` file to dock the ligand into
* `--smiles`: a quoted SMILES str representing the ligand to dock (i.e. "CCC")
* `--hint`: a `.pdb` file with a hint ligand to denote the docking region
* `--hint_radius`: an Angstrom distance from the hint ligand (see above) to consider as the docking region
* `--hint_molinfo`: the resname of the hint ligand used in the hint `.pdb` file
* `--output-dir`: directory to save final docking files (docked ligand and receptor files)

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


