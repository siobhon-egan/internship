#!/bin/bash

mkdir -p ../exout

$usearch -sintax ../out/otus.fa -db ../sintax/rdp_16s_v16.fa -strand plus -tabbedout ../exout/otus_sintax.txt

$usearch -sintax_summary ../exout/otus_sintax.txt -rank p \
  -output ../exout/sintax_summary.txt

cat ../exout/sintax_summary.txt
