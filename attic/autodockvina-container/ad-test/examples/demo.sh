#!/bin/bash

mkdir -p ./out
echo created out directory
/usr/bin/time -p vina --config /examples/config.in --ligand /examples/ZINC00000567.pdbqt --out ./out/ZINC00000567.pdbqt > /dev/null;
score=$(/bin/grep -m 1 "REMARK VINA RESULT:" ./out/ZINC00000567.pdbqt | /usr/bin/awk '{ print $4 }');
echo "./out/ZINC00000567.pdbqt $score" >> ./out/tempResults

cat ./out/tempResults
