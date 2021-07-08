# Description
Calculates fake LogD by outputting the LogP of the given solute and ignoring the solvents using the average of `oechem` and `rdkit`'s LogP calculation

# Setup
1. `cp oe_license.txt SAMPL-league/attic/predict-logd-from-oerdkitlogp`

# Build
1. `cd SAMPL-league/attic/predict-logd-from-oerdkitlogp`
2. `docker build -t oerdlogd .`

# Run
### Options
```
docker run -it --rm oerdlogd --help
Usage: print-LogD [OPTIONS]

  takes in all inputs required for a LogD calculation (solute, solventA and
  solventB) but only calculates the LogP using the average of oechem and
  rdkit logP. Ignores the solventA and solventB inputs

Options:
  --fuzz           Random change logP value by +/- 10%  [default: False]
  --solute TEXT    solute SMILES string
  --solventA TEXT  solventA SMILES string
  --solventB TEXT  solventB SMILES string
  --help           Show this message and exit.
```

### Example Run Commands
`docker ../../ever_given/run.py oerdlogd --solute <solute_SMILES>  --solventA <solventA_SMILES> --solventB <solventB_SMILES>`
* Ex: `docker ../../ever_given/run.py oerdlogd --solute "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  --solventA "O" --solventB "CCCCCCCCO"`
