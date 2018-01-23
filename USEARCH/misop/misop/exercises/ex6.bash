#!/bin/bash

$usearch -sintax_summary ../exout/otus_sintax.txt -rank p \
  -output ../exout/phylum_summary.txt

cat ../exout/phylum_summary.txt
