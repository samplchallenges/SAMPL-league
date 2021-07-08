# Description: 
Outputs LogP of the solute SMILES string using the average of `oechem` and `rdkit`'s LogP calculations

# Setup: 
1. `cp oe_license.txt SAMPL-league/attic/predict-oerdkitlogp`

# Build: 
1. `cd SAMPL-league/attic/predict-oerdkitlogp`
2. `docker build -t oerdlogp .`

# Run:
### Options
```
% docker run -it --rm oerdlogp --help
Usage: print-LogP [OPTIONS]

Options:
  --fuzz             Random change logP value by +/- 10%  [default: False]
  --smiles TEXT      ligand SMILES string
  --output-dir PATH  Output Directory
  --help             Show this message and exit.
```
### Run Command
`python ../../ever_given/run.py oerdlogp --smiles <SMILES_str>`
* `python ../../ever_given/run.py oerdlogp --smiles CCCCc1ccccc1"`
