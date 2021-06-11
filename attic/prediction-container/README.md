# How to build

Copy `oe_license.txt` into this folder (it should look like this when you run the `tree` command)

```
$ tree
.
├── Dockerfile
├── environment.yml
├── oe_license.txt
├── print_logP.py
├── README.md
└── setup.py
```

NOTE: Do not push this container to a public repository or commit the `oe_license.txt` license file (you may push the container to a private repository, we just want to be sure not to leak the file to the public).

Then run `docker build -t calc-logp:0.1 .` to build the container.

# How to use

Run `docker run --rm -it calc-logp:0.1  "CCC"` to calculate the logP of propane.
Run `docker run --rm -it calc-logp:0.1  --help` to see options.
