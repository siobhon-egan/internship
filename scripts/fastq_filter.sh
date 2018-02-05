#!/bin/bash

# This script works with usearch9.2 and usearch10
#
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Quality control and removing dimer seqs
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo ""

		mkdir quality_filered_test

#*****************************************************************************************
for file3 in merging_test/*.fastq
	do
		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Quality control and removing dimer seqs
		echo input is:
		echo ${file3}

	usearch10 -fastq_filter ${file3} -fastaout "quality_filered_test/$(basename "$file3" .fastq).fasta" -fastq_maxee_rate 0.01 -fastq_minlen 150
done
