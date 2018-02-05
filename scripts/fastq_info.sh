#!/bin/bash

raw_seqs="raw_data"
read_summary="1_read_summary"


mkdir $read_summary

 for fq in $raw_data/*R1*.fastq
 do
  usearch10 -fastx_info $fq -output $read_summary/1a_fwd_fastq_info.txt
 done

 for fq in $raw_data/*R2*.fastq
 do
  usearch10 -fastx_info $fq -output $read_summary/1b_rev_fastq_info.txt
 done
