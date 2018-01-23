#!/bin/bash

if [ x$usearch == x ] ; then
	echo Must set \$usearch >> /dev/stderr
	exit 1
fi

version=`$usearch -version | sed "-es/usearch //" | sed "-es/v10.*/v10/"`

if [ x$version != xv10 ] ; then
	echo "usearch version too old, need v10" >> /dev/stderr
	exit 1
fi

out=../out

rm -rf $out
mkdir -p $out

cd $out

# Assemble paired reads, put sample names into read labels
$usearch -fastq_mergepairs ../fq/*_R1_*.fastq -relabel @ -fastqout merged.fq

# Discard reads which probably have errors (quality filtering)
$usearch -fastq_filter merged.fq -fastq_maxee 1.0 -relabel Filt -fastaout filtered.fa

# Find unique sequences and abundances
$usearch -fastx_uniques filtered.fa -sizeout -relabel Uniq -fastaout uniques.fa

# Create 97% OTUs
$usearch -cluster_otus uniques.fa -relabel Otu -otus otus.fa

# Create ZOTUs by denoising (error-correction)
$usearch -unoise3 uniques.fa -zotus zotus.fa

# Create OTU table for 97% OTUs
$usearch -otutab merged.fq -otus otus.fa -strand plus -otutabout otutab.txt

# Create OTU table for ZOTUs
$usearch -otutab merged.fq -zotus zotus.fa -strand plus -otutabout zotutab.txt
