# Setup
Not applicable


# Build
1. `cd SAMPL-league/attic/predict-logd-from-rdkitlogp`
2. `docker build -t logd-calc .`


# Run 
### Options
```
docker run -it --rm logd-calc --help
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
`docker run logd-calc --solute <solute_SMILES> --solventA <solventA_SMILES> --solventB <solventB_SMILES>`
* Ex: `docker run logd-calc --solute "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" --solventA "O" --solventB "CCCCCCCCO"`
