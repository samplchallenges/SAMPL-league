# Description
Calculates fake LogD by outputting the LogP of the given solute and ignoring the solvents using `oechem`'s LogP calculation


# Setup
1. `cp oe_license.txt SAMPL-league/attic/predict-logd-from-oelogp`

# Build
1. `cd SAMPL-league/attic/predict-logd-from-oelogp`
2. `docker build -t oelogd .`


# Run 
### Options
```
docker run -it --rm oelogd --help
Usage: print-LogD [OPTIONS]

  takes in all inputs required for a LogD calculation (solute, solventA and
  solventB) but only calculates the LogP using oechem and ignores the
  solventA  and solventB inputs

Options:
  --fuzz           Randomly change logP value by +/- 10% [default: False]
                   [default: False]

  --solute TEXT    solute SMILES string
  --solventA TEXT  solventA SMILES string
  --solventB TEXT  solventB SMILES string
  --help           Show this message and exit.
  ```
  
### Example Run Commands
`python ../../ever_given/run.py oelogd --solute <solute_SMILES> --solventA <solventA_SMILES> --solventB <solventB_SMILES>`
* Ex: `docker ../../ever_given/run.py oelogd --solute "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" --solventA "O" --solventB "CCCCCCCCO"`

`docker run --rm oelogd <entry_point> --solute <SMILES_str>`
* `docker run --rm oelogd print-LogD --solute "CC(C)CC1=CC=C(C=C1)C(C)C(=O)"`
