#!/bin/bash -x

for dataset in "wiki-vote" "epinions" "wikisigned" "youtube" "pokec" "dbpedia"
do
    for mech in 2 3
    do
        for eps in 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
        do
            ./LDP $dataset $eps 2-9-9 $mech 1.0 1 100
        done
    done
done
