# Description
* This container takes a ligand SMILES string and receptor and hint pdb to specify the binding site, uses `rdkit` to add 3D coordinates, then uses AutoDock Vina to dock the ligand into the specified binding site. 
* This `README.md` details how to setup, build and run the Autodock Vina docker container. 

# Setup:

# Build:
1. `cd adv-base`
2. `docker build -t adv-base .`
3. `cd adv`
4. `docker build -t adv .`


# Run: 
### Options
```
% docker run -it --rm adv --help
Usage: run-autodock [OPTIONS]

  docks the given smiles string into the receptor within the area specified by
  hint and hint-radius

Options:
  --receptor PATH      path of receptor PDB to dock the ligand into
                       [required]
  --smiles TEXT        SMILES str of ligand to be docked. quote to prevent CLI
                       errors "CCC"  [required]
  --hint PATH          path of hint ligand complex for docking region hint
                       [required]
  --hint_molinfo TEXT  residue name of the ligand in the hint complex
                       [required]
  --hint_radius FLOAT  box size of the box to dock into  [required]
  --output-dir TEXT    Output directory for receptor and docked_ligand files
  --debug              prints debug print statements when --debug flag is used
  --help               Show this message and exit.
```

### Example run commands
`python ever_given/run.py adv --file-receptor <receptorPDB> --file-hint <hintPDB> --hint_radius <float> --hint_ligname <str> --smiles <SMILESstr> --output-keys docked_ligand,receptor`
* `python ever_given/run.py adv --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint-ligname E4J --smiles "c1ccc(C(C)C)cc1CNCC(O)(O)[C@@H](NC(=O)[C@@H]2C)C[C@H](C)CCCCCCCCC(=O)N2C" --output-keys docked_ligand,receptor
