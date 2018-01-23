#!/bin/bash

if [ x$usearch == x ] ; then
	echo Must set \$usearch >> /dev/stderr
	exit 1
fi

version=`$usearch -version | sed "-es/usearch //" | sed "-es/v10.*/v10/"`

if [ x$version != xv10 ] ; then
	echo "Wrong usearch version, need v10" >> /dev/stderr
	exit 1
fi

out=../out_mock

rm -rf $out
mkdir -p $out

cd $out

# Merge paired reads
$usearch -fastq_mergepairs ../fq/Mock*_R1_*.fastq -fastqout merged.fq

# Discard reads which probably have errors
$usearch -fastq_filter merged.fq -fastq_maxee 1.0 -relabel Filt -fastqout filtered.fq

# Find unique sequences and abundances
$usearch -fastx_uniques filtered.fq -sizeout -relabel Uniq -fastaout uniques.fa

# Make 97% OTUs
$usearch -cluster_otus uniques.fa -relabel Otu -otus otus.fa

# Error-correction to find all biological sequences (ZOTUs)
$usearch -unoise3 uniques.fa -zotus zotus.fa

# Classify OTUs by aligning to mock reference sequences
$usearch -uparse_ref otus.fa -db ../mockref/HMP_MOCK.v35.fasta -threads 1 \
  -strand both -uparseout otus.uparseref

# Classify ZOTUs by aligning to mock reference sequences
$usearch -uparse_ref zotus.fa -db ../mockref/HMP_MOCK.v35.fasta \
  -threads 1 -strand both -uparseout zotus.uparseref
