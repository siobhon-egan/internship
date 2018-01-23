#!/bin/bash

mkdir -p ../exout

$usearch -otutab_norm ../out/otutab.txt -sample_size 5000 -output ../exout/otutab_5k.txt

$usearch -alpha_div ../exout/otutab_5k.txt -metrics richness,shannon_2 -output ../exout/alpha_5k.txt
