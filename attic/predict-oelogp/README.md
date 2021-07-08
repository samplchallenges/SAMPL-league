# Description:


# Setup:
1. `cp oe_license.txt SAMPL-league/attic/predict-oelogp`


# Build:
1. `cd SAMPL-league/attic/predict-oelogp`
2. `docker build -t oelogp .`


# Run:

### Options
```
% docker run -it --rm oelogp --help
Usage: print-LogP [OPTIONS]

Options:
  --fuzz             Random change logP value by +/- 10%  [default: False]
  --smiles TEXT      ligand SMILES string
  --output-dir PATH  path to output directory
  --help             Show this message and exit.
```

### Run Commands
`python ../../ever_given/run.py oelogp --smiles <SMILES_str>`
* `python ../../ever_given/run.py oelogp --smiles "CCCCc1ccccc1"`
