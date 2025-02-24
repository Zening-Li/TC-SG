#!/bin/bash -x

for dataset in "wiki-vote" "epinions" "wikisigned" "youtube" "pokec" "dbpedia"
do
    for mech in 0 3
    do
        for eps in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5
        do
            ./DPSaveMem $dataset $eps $mech 1.0 1 100
        done
    done
done
