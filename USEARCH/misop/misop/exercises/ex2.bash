#!/bin/bash

mkdir -p ../exout

$usearch -alpha_div ../out/otutab.txt -metrics richness,shannon_2 -output ../exout/alpha.txt
