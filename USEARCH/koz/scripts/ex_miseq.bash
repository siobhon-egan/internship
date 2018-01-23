#!/bin/bash

if [ x$usearch == x ] ; then
	echo Must set \$usearch >> /dev/stderr
	exit 1
fi

rm -rf ../out
mkdir -p ../out
cd ../out

# Merge paired reads
# Add sample name to read label (-relabel option)
# Pool samples together in raw.fq (Linux cat command)
for Sample in Mock Soil Human Mouse
do
	$usearch -fastq_mergepairs ../data/${Sample}*_R1.fq \
	  -fastqout $Sample.merged.fq -relabel $Sample.
	cat $Sample.merged.fq >> all.merged.fq
done

# Strip primers (V4F is 19, V4R is 20)
$usearch -fastx_truncate all.merged.fq -stripleft 19 -stripright 20 \
  -fastqout stripped.fq

# Quality filter
$usearch -fastq_filter stripped.fq -fastq_maxee 1.0 \
  -fastaout filtered.fa -relabel Filt

# Find unique read sequences and abundances
$usearch -fastx_uniques filtered.fa -sizeout -relabel Uniq -fastaout uniques.fa

# Make 97% OTUs and filter chimeras
$usearch -cluster_otus uniques.fa -otus otus.fa -relabel Otu

# Denoise: predict biological sequences and filter chimeras
$usearch -unoise3 uniques.fa -zotus zotus.fa

##################################################
# Downstream analysis of OTU sequences & OTU table
# Can do this for both OTUs and ZOTUs, here do
# just OTUs to keep it simple.
##################################################

# Make OTU table
$usearch -otutab all.merged.fq -otus otus.fa -otutabout otutab_raw.txt

# Normalize to 5k reads / sample
$usearch -otutab_norm otutab_raw.txt -sample_size 5000 -output otutab.txt

# Alpha diversity
$usearch -alpha_div otutab.txt -output alpha.txt

# Make OTU tree
$usearch -cluster_agg otus.fa -treeout otus.tree

# Beta diversity
mkdir beta/
$usearch -beta_div otutab.txt -tree otus.tree -filename_prefix beta/

# Rarefaction
$usearch -alpha_div_rare otutab.txt -output rare.txt

# Predict taxonomy
$usearch -sintax otus.fa -db ../data/rdp_16s_v16.fa -strand both \
  -tabbedout sintax.txt -sintax_cutoff 0.8

# Taxonomy summary reports
$usearch -sintax_summary sintax.txt -otutabin otutab.txt -rank g -output genus_summary.txt
$usearch -sintax_summary sintax.txt -otutabin otutab.txt -rank p -output phylum_summary.txt

# Find OTUs that match mock sequences
$usearch -uparse_ref otus.fa -db ../data/mock_refseqs.fa -strand plus \
  -uparseout uparse_ref.txt -threads 1
