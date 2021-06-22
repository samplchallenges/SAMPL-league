#!/bin/bash
mkdir -p /Users/megosato/SAMPL-league/attic/autodockvina-container/ad-test/out/
/usr/bin/time -p vina --config examples/config.in --ligand examples/ZINC00000567.pdbqt --out /Users/megosato/SAMPL-league/attic/autodockvina-container/ad-test/out/Zinc-out.pdbqt;
