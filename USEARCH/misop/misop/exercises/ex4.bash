#!/bin/bash

mkdir -p ../exout

grep other ../out_mock/otus.uparseref

$usearch -fastx_getseq ../out_mock/otus.fa -label Otu20 -fastaout ../exout/otu20.fa

$usearch -usearch_global ../exout/otu20.fa -db ../out/otus.fa -strand plus -id 0.97 \
  -uc ../exout/otu20_hit.uc

cat ../exout/otu20_hit.uc
