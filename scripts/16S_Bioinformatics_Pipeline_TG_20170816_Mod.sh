#!/bin/bash

# ----------------------------------STEPS FOR 16S 27F-Y/338R V1-2 ANALYSIS FOR 2X250 PAIRED MISEQ READS---------------------------------

# Step 1. Review the length and average error rate of the forward and reverse reads

# Step 2. Merge the paired-end reads.

# Step 3. Trim primers

# Step 4. Quality filtering

# Step 5. Dereplicating sequences

# Step 6. Discarding singletons

# Step 7. Chimera filtering and clustering OTUs

# Step 8. Chimera filtering and clustering ZOTUs

# Step 9. Creating OTU table

# Step 10. Creating ZOTU table

# Step 11. Quality control for OTUs

# Step 12. Assigning taxonomy

# Step 13. Measuring alpha diversity

# Step 14. Measuring beta diversity

# Step 15. Building trees and distance matrices

# --------------------------------------------------------------REFERENCES--------------------------------------------------------------

# SEARCH_16S algorithm
# Edgar, R.C. (2017), SEARCH_16S: A new algorithm for identifying 16S ribosomal RNA genes in contigs and chromosomes.
# http://biorxiv.org/content/early/2017/04/04/124131

# UNBIAS algorithm
# UNBIAS: An attempt to correct abundance bias in 16S sequencing, with limited success.
# http://biorxiv.org/content/early/2017/04/04/124149

# SINAPS algorithm
# SINAPS: Prediction of microbial traits from marker gene sequences. http://biorxiv.org/content/early/2017/04/04/124156

# UNCROSS algorithm
# Edgar, R.C. (2016), UNCROSS: Filtering of high-frequency cross-talk in 16S amplicon reads. doi: http://dx.doi.org/10.1101/088666

# UNOISE algorithm
# Edgar, R.C. (2016), UNOISE2: Improved error-correction for Illumina 16S and ITS amplicon reads.http://dx.doi.org/10.1101/081257

# SINTAX algorithm
# Edgar, R.C. (2016), SINTAX, a simple non-Bayesian taxonomy classifier for 16S and ITS sequences, http://dx.doi.org/10.1101/074161.

# UCHIME2 algorithm
# Edgar, R.C. (2016), UCHIME2: Improved chimera detection for amplicon sequences, http://dx.doi.org/10.1101/074252.

# Expected error filtering and paired read merging
# Edgar, R.C. and Flyvbjerg, H (2014) Error filtering, pair assembly and error correction for next-generation sequencing reads
#   [doi: 10.1093/bioinformatics/btv401].

# UPARSE algorithm
# Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods
# [Pubmed:23955772,  dx.doi.org/10.1038/nmeth.2604].

# USEARCH and UCLUST algorithms
# Edgar,RC (2010) Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461.
# doi: 10.1093/bioinformatics/btq461

# UCHIME algorithm
# Edgar,RC, Haas,BJ, Clemente,JC, Quince,C, Knight,R (2011) UCHIME improves sensitivity and speed of chimera detection,
# Bioinformatics doi: 10.1093/bioinformatics/btr381 [Pubmed 21700674].

# UBLAST algorithm
# Edgar, RC (unpublished) http://drive5.com/usearch


# ----------------------------------------------------------------FOLDERS----------------------------------------------------------------

usearch10="usearch10.0.240_i86osx64"
usearch9="usearch10.0.240_i86osx64"
raw_seqs="raw_data"
read_summary="1_read_summary"
merged_seqs="2a_merged_seqs"
fwd_not_merged="2b_fwd_not_merged"
rev_not_merged="2c_rev_not_merged"
trimmed_seqs="3_trimmed_output"
fasta_seqs_filtered="4a_fasta_seqs_filtered"
fasta_seqs_discarded="4b_fasta_seqs_discarded"
unique_seqs="5_unique_seqs"
no_singletons="6_no_singletons"
OTUs="7_OTUs"
ZOTUs="8_ZOTUs"
otu_table="9_otu_table"
zotu_table="10_zotu_table"
quality_control_otus="11a_quality_control_otus"
quality_control_zotus="11b_quality_control_zotus"
assigned_taxonomy_otus="12a_assigned_taxonomy_otus"
assigned_taxonomy_zotus="12b_assigned_taxonomy_zotus"
alpha_diversity="13_alpha_diversity"
beta_diversity_otus="14a_beta_diversity_OTUs"
beta_diversity_zotus="14b_beta_diversity_ZOTUs"
trees_otus="15a_Trees_OTUs"
trees_zotus="15b_Trees_ZOTUs"

# --------------------------------------16S 27F-Y/338R V1-2 ANALYSIS FOR 2X250 PAIRED MISEQ READS--------------------------------------

echo ""
figlet Starting 16S Analysis
echo ""

# -------------------------------------------------------------READ SUMMARY-------------------------------------------------------------

# Read summary

# The fastx_info command is very useful for getting a quick summary of the length distribution and quality of the reads in a FASTQ file.

# Check all of the FASTQ files in case some of them are different for some reason.

# Use grep to get a summary of one thing at a time (expected error (EE) distribution, length distribution, number of reads...) to check
# that all of the files are similar, otherwise some may need special treatment. E.g.,
# cd ../fastq_info
# grep "^EE" *

# If you have paired reads, review the R1 and R2 files separately because they are usually different, e.g. the R2s tend to have lower
# quality.

# Use the fastq_eestats2 command to review whether there are low-quality tails that should be truncated. If you have reads that vary in
# length, e.g. from 454 pyrosequencing, fastq_eestats2 will help you to decide which length to choose for trimming.

# Get the primer sequences (the segments of the PCR oligos that bind to your gene, e.g. 16S primers), and confirm that at least the
# forward primer is present in the reads. If you have paired reads, then the R1s usually start with the forward primer and the R2s
# usually start with the reverse primer. You can do this using the search_oligodb command, e.g.

# usearch10.0.240_i86osx32 -search_oligodb otus.fa -db primers.fa -strand both \
#  -userout primer_hits.txt -userfields query+qlo+qhi+qstrand

# If you don't find the primers then you should make sure you understand why. Maybe the library was prepared using an unusual strategy
# such as nested PCR, which may need some special processing in the pipeline, or maybe you have the wrong primer sequences.

# If the FASTQ files are large, you can use the fastx_subsample command to get a smaller subset of reads for quick testing.

echo ""
figlet -f bubble Step 1
figlet -f small Read Summary
echo ""

echo ""
figlet -f digital CREATING READ SUMMARY DIRECTORY
echo ""

mkdir $read_summary

echo ""
figlet -f digital FORWARD READ SUMMARY
echo ""
 for fq in $raw_seqs/*R1*.fastq
 do
  $usearch10 -fastx_info $fq -output $read_summary/1a_fwd_fastq_info.txt
 done
echo ""
figlet -f digital REVERSE READ SUMMARY
echo ""
 for fq in $raw_seqs/*R2*.fastq
 do
  $usearch10 -fastx_info $fq -output $read_summary/1b_rev_fastq_info.txt
 done
echo ""
figlet -f digital OUTPUT FILES:
echo $read_summary
echo 1a_fwd_fastq_info.txt
echo 1b_rev_fastq_info.txt
echo ""

# -------------------------------------------------------MERGING PAIRED-END READS-------------------------------------------------------

echo ""
figlet -f bubble Step 2
figlet -f small Merging paired-end reads
echo ""

# forward reads: -fastq_mergepairs *R1*.fastq

# reverse reads: -reverse *R2*.fastq

# output files: -fastqout merged.fq

# If the -reverse option is omitted, the reverse FASTQ filename is constructed by replacing R1 with R2.

# You can use shell wildcards (* and ?) to give a pattern that matches the FASTQ files you want to merge.

# If multiple samples are combined into a single file as shown in some of the above examples, then you lose track of which read came from
# which sample. This is addressed by adding a sample identifier to each read label. The simplest method is to use the -sample option, e.g.
# usearch -fastq_mergepairs SampleA_R1.fastq -fastqout merged.fq -sample SampleA. The string sample=SampleA; will be added at the end of
# the read label.

# FASTQ filenames are often based on the sample identifier, e.g. SampleA_R1.fastq. If you specify  -relabel @ then fastq_mergepairs gets
# the sample identifier from the FASTQ file name by truncating at the first underscore (_) or period (.). A period and the read number is
# added after the sample identifier to make the new read label, which replaces the original label. This differs from the -sample option,
# which adds the sample= annotation at the end of the label. The usearch_global command understands both of these methods for putting
# sample identifiers into read labels.. usearch -fastq_mergepairs SampleA_R1.fastq -fastqout merged.fq  -relabel @

#  -fastqout_notmerged_fwd  FASTQ filename for forward reads which were not merged
#  -fastqout_notmerged_rev  FASTQ filename for reverse reads which were not merged

# Reports
#  -report   Filename for summary report. See Reviewing a fastq_mergepairs report to check for problems.
#  -tabbedout  Tabbed text file containing detailed information about merging process for each pair including reason for discarding.
#  -alnout  Human-readable alignments. Useful for trouble-shooting.
#
# Merged read labels
#  -relabel  Prefix string for output labels. The read number 1, 2, 3... is appended after the prefix.
#  -relabel @ Relabel using prefix string constructed from FASTQ filename, this will be understood as the sample identifier.
#  -sample  xxx Append sample identifier to read label using sample=xxx; format. This is an alternative method for adding sample ids.
#  -fastq_eeout  Add ee=xxx; annotation with the number of expected errors in the merged read.
#  -label_suffix  Suffix to append to merged read label. Can be used e.g. to add sample=xxx; type of sample identifier annotations.
#
# Filtering
#  -fastq_maxdiffs  Maximum number of mismatches in the alignment. Default 5. Consider increasing if you have long overlaps.
#  -fastq_pctid  Minimum %id of alignment. Default 90. Consider decreasing if you have long overlaps.
#  -fastq_nostagger  Discard staggered pairs. Default is to trim overhangs (non-biological sequence).
#  -fastq_merge_maxee  Maximum expected errors in the merged read. Not recommended for OTU analysis.
#  -fastq_minmergelen  Minimum length for the merged sequence. See Filtering artifacts by setting a merge length range.
#  -fastq_maxmergelen  Maximum length for the merged sequence.
#  -fastq_minqual  Discard merged read if any merged Q score is less than the given value. (No minimum by default).
#  -fastq_minovlen  Discard pair if alignment is shorter than given value. Default 16.
#
# Pre-processing of reads before alignment
# -fastq_trunctail  Truncate reads at the first Q score with <= this value. Default 2.
# -fastq_minlen  Discard pair if either read is shorter than this, after truncating by -fastq_trunctail if applicable. Default 64.
#
# Multi-threading
# -threads  Specifies the number of threads. Default 10, or the number of CPU cores, which ever is less.

echo ""
figlet -f digital CREATING MERGED AND UNMERGED SEQUENCE DIRECTORIES
echo ""

mkdir $merged_seqs
mkdir $fwd_not_merged
mkdir $rev_not_merged

figlet -f digital MERGING PAIRED-END READS

usearch10 -fastq_mergepairs raw_data/*R1*.fastq -reverse raw_data/*R2*.fastq -relabel @ -fastq_maxdiffs 15 -fastq_minovlen 20 \
-fastq_minmergelen 30 -fastq_maxmergelen 500 -fastqout merged_seqs/merged.fq -fastqout_notmerged_fwd fwd_not_merged/not_merged_fwd.fq \
-fastqout_notmerged_rev rev_not_merged/not_merged_rev.fq -report merged_seqs/2d_merging_seqs_report.txt -tabbedout merged_seqs/2e_tabbedout.txt \
-alnout merged_seqs/2f_aln.txt

echo ""
figlet -f digital OUTPUT FILES:
echo $merged_seqs
echo $fwd_not_merged
echo $rev_not_merged
echo 2d_merging_seqs_report.txt
echo 2e_tabbedout.txt
echo 2f_aln.txt
echo ""

# ----------------------------------------------------------TRIMMING PRIMERS----------------------------------------------------------

# fastx_truncate command
# Commands > FASTA/Q files, Reads

# If the reads start with a binding sequence for a primer of length L, then a simple way to strip is to delete the first L
# letters of the reads using the fastx_truncate command.

# Truncates and / or pads sequences in a FASTA or FASTQ.
# If the -trunclen n option is specified, sequences are truncated at length n. Shorter sequences are discarded.
# If the -padlen option is specified, padding occurs before truncating.

# Forward primer: Bact_16S_27F-Y_v1-v2_primer - AGAGTTTGATCCTGGCTYAG = 20 bases
# Reverse primer: Bact_16S_338R_v1-v2_primer - TGCTGCCTCCCGTAGGAGT = 19 bases

# Output is written to -fastaout (FASTA) and/or -fastqout (FASTQ).
# The -stripleft n option removes n letters at the start of a sequence.
# The -stripright n option removes n letters at the end of a sequence.
# The -padlen n option appends Ns at the end of the sequence as needed to increase the length to n. If the length is already >= n,
# then the sequence is unchanged. (Requires v8.0.1611 or later).
# The -padq c option specifies the ASCII character to use for padding the quality scores in a FASTQ file if -padlen is specified.
# Default is I (upper-case letter i), which means Phred score Q=40 with most current FASTQ files (see FASTQ file parameters for
# discussion).
# The -relabel xxx option specifies that new sequence labels should be generated. The new labels are xxx followed by 1, 2, 3...

echo ""
figlet -f bubble Step 3
figlet -f small Trimming primers
echo ""

echo ""
figlet -f digital CREATING PRIMER-TRIMMED SEQUENCE DIRECTORY
echo ""

mkdir $trimmed_seqs

$usearch10 -fastx_truncate $merged_seqs/merged.fq -stripleft 20 -stripright 19 -fastqout $trimmed_seqs/trimmed.fq

echo ""
figlet -f digital OUTPUT FILES:
echo $trimmed_seqs
echo ""

# ----------------------------------------------------------QUALITY FILTERING----------------------------------------------------------
# fastq_filter command
# Commands > Reads

# Performs quality filtering and / or conversion of a FASTQ file to FASTA format.
# See also
#   Paper describing expected error filtering and paired read merging (Edgar & Flyvbjerg, 2015).
#   Paired read assembler and quality filtering benchmark results
#   Read quality filtering
#   Expected errors
#   FASTQ format options
#   Quality scores
#   Choosing FASTQ filter parameters
#   Strategies for dealing with low-quality reverse reads (R2s)

# The fastx_learn command is useful for checking the error rate after expected error quality filtering, which assumes that the Q scores
# are accurate. It does not use Q scores so gives an independent check.

# Output options

# Option	 	Description
# -fastqout filename	 	FASTQ output file. You can use both ‑fastqout and ‑fastaout.
# -fastaout filename	 	FASTA output file. You can use both ‑fastqout and ‑fastaout.
# -fastqout_discarded filename	 	FASTQ output file for discarded reads. You can use both ‑fastqout_discarded and ‑fastaout_discarded.
# -fastaout_discarded filename	 	FASTA output file for discarded reads. You can use both ‑fastqout_discarded and ‑fastaout_discarded.
# -relabel prefix	 	Generate new labels for the output sequences. They will be labeled prefix1, prefix2 and so on. For example, if
# you use -relabel SampleA. then the labels will be SampleA.1, SampleA.2 etc.
# The special value @ indicates that the string should be constructed from the file name by truncating the file name at the first
# underscore or period and appending a period. With a typical Illumina FASTQ file name, this gives the sample name. So, for example, if
# the FASTQ file name is Mock_S188_L001_R1_001.fastq, then the string is Mock and the output labels will be Mock.1, Mock.2 etc.
# -fastq_eeout	 	Append the expected number of errors according to the Q scores to the label in the format "ee=xx;". Expected errors
# are calculated after truncation, if applicable.
# -sample string	 	Append sample=string; to the read label.
# Filtering options

# Option	 	Description
# ‑fastq_truncqual N	 	Truncate the read at the first position having quality score <= N, so that all remaining Q scores are >N.
# -fastq_maxee E	 	Discard reads with > E total expected errors for all bases in the read after any truncation options have been
# applied.
# ‑fastq_trunclen L	 	Truncate sequences at the L'th base. If the sequence is shorter than L, discard.
# -fastq_minlen L	 	Discard sequences with < L letters.
# -fastq_stripleft N	 	Delete the first N bases in the read.
# -fastq_maxee_rate E	 	Discard reads with > E expected errors per base.Calculated after any truncation options have been applied.
# For example, with the fastq_maxee_rate option set to 0.01, then a read of length 100 will be discarded if the expected errors is >1,
# and a read of length 1,000 will be discarded if the expected errors is >10.
# -fastq_maxns k	 	Discard if there are >k Ns in the read.

# Examples

# "Raw" conversion of Sanger FASTQ to FASTA with no filtering:

#   usearch -fastq_filter reads.fastq -fastaout reads.fasta -fastq_ascii 64

# Truncate to length 150, discard if expected errors > 0.5, and convert to FASTA:
#
#   usearch -fastq_filter reads.fastq -fastq_trunclen 150 -fastq_maxee 0.5 \
#     -fastaout reads.fasta

# Truncate a read at length 100 and then discard if it contains a Q<15, output to new FASTQ file:

#   usearch -fastq_filter reads.fastq -fastq_minlen 100 -fastq_truncqual 15 \
#     -fastqout filtered_reads.fastq

# Quality (Phred) scores

# See also
#   FASTQ files
#   Average Q is a bad idea!
#   Expected errors
#   Quality filtering
# The quality score of a base, also known as a Phred or Q score, is an integer value representing the estimated probability of an error,
# i.e. that the base is incorrect. If P is the error probability, then:
#
#   P = 10–Q/10
#
#   Q = –10 log10(P)

# Q scores are often represented as ASCII characters. The rule for converting an ASCII character to an integer varies, see FASTQ options
# for details. Tables converting between integer Q scores, ASCII characters and error probabilities are shown in the table below
# ASCII_BASE 33, which is now almost universally used, and ASCII_BASE 64 which is used in some older Illumina data.

# What kind of error?
# There is an important difference between Q scores in reads from 454 and Illumina. In effect, 454 ignores the possibility of
# substitution errors and Illumina ignores indels. With 454, the Q score is the estimated probability that the length of the homopolymer
# is wrong, and with Illumina the Q score is the probability that the base call is incorrect. In the case of Illumina, this is
# reasonable because indel errors are very rare. But with 454, substitution errors are quite common, occurring with comparable
# frequency to homopolymer errors. This means that 454 Q scores are not as informative as Illumina Q scores, but are still useful in
# practice. See quality filtering for further discussion.

# Small Q scores
# Note that a Q score of 3 means P=0.5, meaning that there is a 50% chance the base is wrong, and lower values represent even higher
# probabilities of error. Q=0 means P=1, i.e. that the base call is certainly wrong, so this is rarely used, though might be appropriate
# for an undetermined base (often represented as 'N'). I have never seen a FASTQ file with Q=0, but since the format is not standardized
# I can't be sure. The lowest value usually found in practice is Q=2 (P=0.63), which means the base call is more likely to be wrong than
# correct.

# Recognizing the format
# The fastx_info and fastq_chars commands can be used to determin the format. The most important parameter is ASCII_BASE, which as far
# as I know is always 33 or 64. With a typical range from Q2 to Q40, this gives a range of ASCII values from 35 to 73 with ASCII_BASE=33
# and from 66 to 104 with ASCII_BASE=64. These ranges overlap from ASCII 66 to 73. Also, values >Q40 may be produced by some machine
# software and by some post-processing software such as paired read assemblers. So if we see ASCII values >73 that doesn't necessarily
# mean that we have ASCII_BASE=64, these could be high quality scores with ASCII_BASE=33. The only sure way to distinguish for sure is
# if we see ASCII values < 64, in which case we know ASCII_BASE=33. A quick way to check visually is to look for # and $, which means
# ASCII_BASE=33 or lower-case letters which probably implies ASCII_BASE=64.

# Also see this page for selectiong and E score: http://www.drive5.com/usearch/manual/exp_errs.html

echo ""
figlet -f bubble Step 4
figlet -f small Quality filtering
echo ""

echo ""
figlet -f digital CREATING FILTERED SEQUENCE DIRECTORIES
echo ""

mkdir $fasta_seqs_filtered
mkdir $fasta_seqs_discarded

$usearch10 -fastq_filter $trimmed_seqs/*.fq -fastq_maxee 0.25 -fastq_minlen 250 -fastaout $fasta_seqs_filtered/filtered.fa \
-fastaout_discarded $fasta_seqs_discarded/discarded.fa

echo ""
figlet -f digital OUTPUT FILES:
echo $fasta_seqs_filtered
echo $fasta_seqs_discarded
echo ""

# -------------------------------------------------------DEREPLICATING SEQUENCES-------------------------------------------------------

# fastx_uniques command
# Commands > Clustering, FASTA/Q files, Reads

# See also
#   search_exact command

# Find the set of unique sequences in an input file, also called dereplication. Input is a FASTA or FASTQ file. Sequences are compared
# letter-by-letter and must be identical over the full length of both sequences (substrings do not match). Case is ignored, so an
# upper-case letter matches a lower-case letter.

# All 26 letters of the English alphabet are treated in the same way, so there is no concept of a biological alphabet or of wildcard
# matching (unless strand -both is used).
# Multithreading is supported.

# The -fastaout option specifies a FASTA output file for the unique sequences. Sequences are sorted by decreasing abundance.

# The -fastqout option specifies a FASTQ output file for the unique sequences. Sequences are sorted by decreasing abundance.

# The -tabbedout option specifies an output file in tabbed text format. The fields are: 1. input label, 2. output label (this is the
# input label of the first occurrence of the sequence, or the new label assigned to it if the -relabel option is used), 3. cluster
# number (zero-based, so 0 is the first unique sequence found, 1 is the second etc.), 4. member number in the cluster (zero -based),
# 5. input label of the first occurrence of the sequence (only if -relabel is specified).

# The -uc output file is supported, but not other standard output files.

# The -sizeout option specifies that size annotations should be added to the output sequence labels.

# The -relabel option specifies a string that is used to re-label the dereplicated sequences. An integer is appended to the label.
# E.g., -relabel Uniq will generate sequences labels Uniq1, Uniq2 ... etc. By default, the label of the first occurrence of the
# sequence is used.

# The -minuniquesize option sets a minimum abundance. Unique sequences with a lower abundance are discarded. Default is 1, which means
# that all unique sequences are output.

# The -topn N option specifies that only the first N sequences in order of decreasing abundance will be written to the output file.

# Reverse-complemented matching for nucleotide sequences is supported by specifying -strand both.

# Example

# usearch -fastx_uniques input.fasta -fastaout uniques.fasta -sizeout -relabel Uniq

echo ""
figlet -f bubble Step 5
figlet -f small Dereplicating sequences
echo ""

echo ""
figlet -f digital CREATING FILTERED SEQUENCE DIRECTORIES
echo ""

mkdir $unique_seqs

$usearch10 -fastx_uniques $fasta_seqs_filtered/*.fa -fastaout $unique_seqs/unique.fa -sizeout -relabel Uniq

echo ""
figlet -f digital OUTPUT FILES:
echo $fasta_seqs_filtered
echo $fasta_seqs_discarded
echo ""

# -------------------------------------------------------DISCARDING SINGLETONS-------------------------------------------------------

echo ""
figlet -f bubble Step 6
figlet -f small Discarding singletons
echo ""

echo ""
figlet -f digital CREATING FILTERED SEQUENCE DIRECTORIES
echo ""

mkdir $no_singletons

$usearch10 -sortbysize $unique_seqs/unique.fa -fastaout $no_singletons/no_singletons.fa -minsize 2

echo ""
figlet -f digital OUTPUT FILES:
echo $no_singletons
echo ""

# -------------------------------------------------CHIMERA FILTERING AND CLUSTERING OTUS-------------------------------------------------

# UCHIME algorithm

# See also
#   Chimera formation
#   UCHIME order dependency
#   Reproducibility of results
#   Comparison of UCHIME and DECIPHER
#   Downloads
# UCHIME and UCHIME2 are algorithms for detecting chimeric sequences. See the UCHIME and UCHIME2 papers for details. Note that the
# accuracy results reported in the original UCHIME paper are highly over-optimistic due to unrelaistic benchmark design. This is
# explained in the UCHIME2 paper.

# Stand-alone UCHIME[2] should not be used in OTU clustering or denoising pipelines because the cluster_otus and unoise3 commands
# have more effective chimera filters integrated into the algorithms. You will probably get much higher error rates, both false
# positives and false negatives, if you use UCHIME[2] as a stand-alone filtering step.

# cluster_otus command
# Commands > OTU analysis, Clustering, Chimeras
# See also
#   OTU / denoising pipeline
#   Defining and interpeting OTUs
#   OTU benchmark results
#   Making an OTU table (otutab command)
#   Should I use UPARSE or UNOISE?

# The cluster_otus command performs 97% OTU clustering using the UPARSE-OTU algorithm.

# Chimeras are filtered by this command. This chimera filtering is much better than using UCHIME so I do not recommend using
# reference-based chimera filtering as a post-processing step, except as a manual check, because false positives are common.

# For most purposes, I consider 97% OTU clustering obsolete. It is better to use the unoise command to recover the full set of biological
#  sequences in the reads. These are also valid OTUs; I call them "ZOTUs" for zero-radius OTUs, to emphasize this. See defining and
# interpreting OTUs and the UNOISE paper for further discussions.

# Input to cluster_otus is a FASTA file containing quality filtered, globally trimmed and dereplicated reads from a marker gene amplicon
# sequencing experiment, e.g. 16S or ITS. It is generally recommended that singleton reads should be discarded.

# See OTU / denoising pipeline for discussion of how to prepare reads before clustering. It is strongly recommended that you follow the
# pipeline recommentations, otherwise the accuracy of the OTUs will probably be compromised.

# Input sequence labels must have size annotations giving the abundance of the unique sequence. Size annotations are generated by the
# -sizeout option of clustering commands; typically fastx_uniques is used.

# The -minsize option  can be used to specify a minimum abundance. Default value is 2, which discards singleton unique sequences.

# The identity threshold is fixed at 97%. See defining and intepreting OTUs for discussion. See UPARSE OTU radius for making OTUs at
# different identities.

# The -otus option specifies a FASTA output file for the OTU representative sequences. By default, OTUs labels are taken from the input
# file, with size annotations stripped.

# The -relabel option specifies a string that is used to re-label OTUs. If -relabel xxx is specified, then the labels are xxx followed
# by 1, 2 ... up to the number of OTUs. OTU identifiers in the labels is required for making an OTU table with the otutab command.

# The -uparseout option specifies a tabbed text output file documenting how the input sequences were classified.

# The -uparsealnout option species a text file containing a human-readable alignment of each query sequence to its UPARSE-REF model.

# If you want size=nnn; annotations in the OTU or ZOTU sequence labels, see adding sizes.

# Parsimony score options are supported.

# Alignment parameters and heuristics are supported.

# Example

# usearch -cluster_otus uniques.fa -otus otus.fa -uparseout uparse.txt -relabel Otu

# FAQ: Should you use UPARSE or UNOISE?
# There are two different ways to make OTUs: 97% clustering and denoising.

# The UPARSE algorithm makes 97% OTUs.

# The UNOISE algorithm does denoising, i.e. error-correction.

# If UPARSE works perfectly, it will give you a subset of the correct biological sequences in your reads such that no two sequences
# are >97% identical to each other. It is implemented in the cluster_otus command.

# If UNOISE works perfectly, it will give you all the correct biological sequences in the reads. It is implemented in the unoise3
# command.

# (Of course, no algorithm is perfect so you should expect some mistakes).

# The UPARSE pipeline and UNOISE pipeline are very similar, the main difference is whether you run cluster_otus or unoise3 as the
# clustering step.

# Once you have made an OTU table, you can proceed with diversity analysis etc. in the same way, regardless of whether you used
# UPARSE or UNOISE.

# Which should you choose? I suggest you try both.

# Pros and cons
# Almost all published papers use 97% clustering, so this will be easier to explain to your PI and to referees. The main disadvantage
# of 97% clustering is that you discard some correct biological sequences that are present in your reads. If these represent strains
# or species with a different phenotype, then you lose relevant information and the corresponding reads will be lumped together into
# one OTU that contains multiple phenotypes.

# The main advantage of denoising is that you get better resolution by keeping all the biological sequences. The main disadvantage of
# denoising is that species often have variations between individuals and paralogs that are not 100% identical. If you have
# intra-species variations in the region that you sequenced, then you will get two or more OTUs for one species. For most purposes,
# this really doesn't matter -- it might even be better if this enables you to detect strains with different phenotypes -- so if I have
# to recommend one method, then I would recommend denoising.

echo ""
figlet -f bubble Step 7
figlet -f small Filtering chimeras and clustering OTUs
echo ""

echo ""
figlet -f digital CREATING OTU DIRECTORY
echo ""

mkdir $OTUs

$usearch10 -cluster_otus $no_singletons/no_singletons.fa -otus $OTUs/otus.fa -uparseout $OTUs/7_uparse.txt -relabel OTU

echo ""
figlet -f digital OUTPUT FILES:
echo $OTUs
echo 7_uparse.txt
echo ""

# ------------------------------------------------CHIMERA FILTERING AND CLUSTERING ZOTUS------------------------------------------------

# unoise3 command
# Commands > OTU analysis, Reads, Chimeras
# See also
#   OTU / denoising analysis pipeline
#   UNOISE paper
#   Should I use UPARSE or UNOISE?

# Uses the UNOISE algorithm to perform denoising (error-correction) of amplicon reads.

# Errors are corrected as follows:
#   - Reads with sequencing error are identified and corrected.
#   - Chimeras are removed.

# Input is a set of quality-filtered unique read sequences with size=nnn; abundance annotations. See OTU / denoising pipeline for
# details of how reads should be pre-processed and how other types of errors and artifacts can be removed.

# The algorithm is designed for Illumina reads, it does not work as well on 454, Ion Torrent or PacBio reads.

# Predicted correct biological sequences are written to the -zotus file in FASTA format. Labels are formatted as Zotunnn;Uniqlabel;
# where nnn is 1, 2, 3... and Uniqlabel is the label from the input file (truncated at the first semi-colon, to strip any annotations).

# Predicted correct amplicon sequences are written to the -ampout fle in FASTA format. These include chimeras, so this output file is
# not generally needed in a production pipeline. Labels are formatted as Ampnnn;uniq=Uniqlabel;uniqsize=u;size=s; where nnn is 1, 2,
# 3..., Uniqulabel is the label in the input file, truncated at the first semi-colon, u is the size= annotation from the input file and
# s is the total size of reads derived from this amplicon.

# The -chimeras option species a FASTA file for amplicons which are predicted to be chimeric.

# An OTU table can be generated using the otutab command. See OTU / denoising pipeline.

# The -minsize option specifies the minimum abundance (size= annotation). Default is 8. Input sequences with lower abundances are
# discarded. Most of the low-abundance sequences are usually noisy and are be mapped to a ZOTU by the otutab command. For higher
# sensivity, reducing minsize to 4 is reasonable, especially if samples are denoised indivudually rather pooling all samples together,
# as I would usually recommend. With smaller minsize, there tends to be more errors in the predicted low-abundance biological sequences.

# The -tabbedout option specifies a tabbed text filename which reports the processing done for each sequence, e.g. if it is classified
# as noisy or chimeric.

# The -unoise_alpha option specifies the alpha parameter (see UNOISE2 paper for definition). Default is 2.0.

# If you want size=nnn; annotations in the OTU or ZOTU sequence labels, see adding sizes.

# Example

# usearch -unoise3 uniques.fa -zotus zotus.fa -tabbedout unoise3.txt

echo ""
figlet -f bubble Step 8
figlet -f small Filtering chimeras and clustering ZOTUs
echo ""

echo ""
figlet -f digital CREATING ZOTU DIRECTORY
echo ""

mkdir $ZOTUs

$usearch10 -unoise3 $no_singletons/no_singletons.fa -zotus $ZOTUs/zotus.fa -tabbedout $ZOTUs/8_unoise3.txt

echo ""
figlet -f digital OUTPUT FILES:
echo $ZOTUs
echo 8_unoise3.txt
echo ""

# --------------------------------------------------------CONSTRUCTING OTU TABLE--------------------------------------------------------

# Creating an OTU table
# See also
#   OTU / denoising pipeline
#   Defining and interpreting OTUs
#   Interpreting counts and frequencies in a OTU table
#   Mapping reads to OTUs

# In my approach, OTUs are sequences selected from the reads. The goal is to identify a set of correct biological sequences; see
# defining and interpreting OTUs for discussion.

# An OTU table is a matrix that gives the number of reads per sample per OTU. One entry in the table is usually a number of reads,
# also called a "count", or a frequency in the range 0.0 to 1.0.  Counts are calculated by mapping reads to OTUs. Note that read
# abundance has very low correlation with species abundance, so the biological interpretation of these counts is unclear.

# An OTU table is created by the otutab command. Input should be reads before quality filtering and before discarding low-abundance
# unique sequences, e.g. singletons. This will improve sensitivity, because many low-quality and low-abundance reads will map
# successfully to OTUs. Reads must have sample indenfiers in the labels. The search database is the OTU (or ZOTU) sequences, which must
# have valid OTU identifiers in the labels.

# I suggest making OTU tables for both 97% OTUs and denoised sequences (ZOTUs).

# Normalizing the table
# After generating the table, you should use the otutab_norm command to normalize all samples to the same number of reads.

# Examples

# usearch -otutab reads.fq -otus otus.fa -otutabout otutab.txt -mapout map.txt

# usearch -otutab reads.fq -zotus zotus.fa -otutabout zotutab.txt -mapout zmap.txt

# usearch -otutab_norm otutab.txt -sample_size 5000 -output otutab_norm.txt

# usearch -otutab_norm zotutab.txt -sample_size 5000 -output zotutab_norm.txt

# otutab command
# Commands > OTU analysis
# See also
#   Defining and interpreting OTUs
#   OTU clustering
#   UPARSE pipeline
#   UNOISE pipeline
#   OTU commands

# The otutab command generates an OTU table by mapping reads to OTUs.

# OTU table output
# See OTU table output options.

# Normalizing the table
# After generating the table, you should use the otutab_norm command to normalize all samples to the same number of reads.

# Query dataset
# The query file can be in FASTQ or FASTA format. Every query sequence must be labeled with a sample identifier. The
# fastx_get_sample_names command can be used to check that your sample names are formatted correctly.

# fastx_get_sample_names command
# Reports the sample identifiers found in seqence labels. The -output option specifies a text file to contain one line per identifier.
# This is useful as a check that your sample identifiers are formatted correctly.

# Example

# usearch -fastx_get_sample_names reads.fa -output samples.txt

# Query sequences are typically raw reads, i.e. reads after paired read merging, if applicable, but before quality filtering.
# Low-quality reads and singletons can often be mapped successfully to an OTU, so including them accounts for a larger fraction of the
# reads (see OTU coverage). The fastx_uniques_persample command can be used to find the unique sequences and abundances for all samples.
# This compresses the input data and makes the otutab command somewhat faster but probably not as much as you might expect
# (typically, the compression is only ~2x).

# OTU database
# The search database is either a set of OTU sequences or "ZOTU" sequences, i.e. denoised sequences. Each query sequence is mapped to
# the closest database sequence. Ties are broken systematically by picking the first in database file order. A udb database can be used.
# Database sequences must be labeled with OTU identifiers. The database file is specified by the -otus or -zotus option. Use -zotus if
# the OTUs are denoised, -otus otherwise.

# OTU sizes
# See adding sizes to OTU labels.

# Identity threshold for mapping
# The -id option sets the minimum fractional identity. Default is 0.97, corresponding to 97% identity. Denoised OTUs also use a
# 97% identity threshold by default to allow for sequencing and PCR error. See defining and interpreting OTUs and mapping reads to
# OTUs for discussion.

# By default, reads are assumed to be on the same strand as the OTU sequences. You can use -strand both to search both strands.

# The -notmatched option specifies a FASTA filename for sequences which are not assigned to an OTU.

# The -notmatchedfq option specifies a FASTQ file for unassigned sequences (input must be FASTQ).

# Multithreading and standard output files are supported.

# Example

# usearch -otutab reads.fq -otus otus.udb -otutabout otutab.txt -biomout otutab.json \
#   -mapout map.txt -notmatched unmapped.fa -dbmatched otus_with_sizes.fa -sizeout

# OTU table output files
# See also
#   OTU / denoising pipeline
#   Generating an OTU table

# Commands which support generating OTU tables understand the following options. The most commonly used commands are otutab and
# closed_ref.

# -mapout filename
#     Tab-separated text file which maps reads to OTUs.
#     There is one read per line with two fields: the read label and the OTU identifier.

# -otutabout filename
#     QIIME classic tabbed text format.

# -biomout filename
#     BIOM v1.0 format (JSON). The biom utility can be used to convert to BIOM v2.1 format (HDF5).

# -mothur_shared_out filename
#     Mothur "shared" file.

echo ""
figlet -f bubble Step 9
figlet -f small Creating OTU table
echo ""

echo ""
figlet -f digital CREATING OTU TABLE DIRECTORY
echo ""

mkdir $otu_table
mkdir $zotu_table

$usearch10 -otutab $trimmed_seqs/trimmed.fq -otus $OTUs/otus.fa -otutabout $otu_table/9_otutab.txt -biomout $otu_table/9_otutab.json \
-mapout $otu_table/9_map.txt -notmatchedfq $otu_table/9_unmapped.fq -matched $otu_table/9_otus_with_sizes.fa -sizeout

echo ""
figlet -f digital OUTPUT FILES:
echo 9_otutab.txt
echo 9_otutab.json
echo 9_map.txt
echo 9_unmapped.fa
echo 9_otus_with_sizes.fa
echo ""

# --------------------------------------------------------CONSTRUCTING ZOTU TABLE--------------------------------------------------------

echo ""
figlet -f bubble Step 10
figlet -f small Creating ZOTU table
echo ""

echo ""
figlet -f digital CREATING ZOTU TABLE DIRECTORY
echo ""

$usearch10 -otutab $trimmed_seqs/trimmed.fq -zotus $ZOTUs/zotus.fa -otutabout $zotu_table/10_zotutab.txt \
-biomout $zotu_table/10_zotutab.json -mapout $zotu_table/10_map.txt -notmatchedfq $zotu_table/10_unmapped.fq \
-dbmatched $zotu_table/10_zotus_with_sizes.fa -sizeout -keep_annots

echo ""
figlet -f digital OUTPUT FILES:
echo 10_zotutab.txt
echo 10_zotutab.json
echo 10_map.txt
echo 10_unmapped.fa
echo 10_zotus_with_sizes.fa
echo ""

# ------------------------------------------------------QUALITY CONTROL FOR OTUS-------------------------------------------------------

# OTU coverage
# See also
#   OTU / denoising analysis pipeline
#   Quality control for OTUs
#   Comparing two sets of OTUs for the same reads

# Explaining your data
# Measuring coverage attempts to answer the question "what fraction of the reads is explained by my OTUs?". This gives a sense of whether
# the OTUs have good coverage. For example, if 95% of your reads can be explained by the OTUs, then you might feel that the coverage is
# good. The 5% of reads which are not explained probably contain a mix of rare biological sequences and artifacts such as sequencer
# errors and chimeras, but with low-abundance sequences these cannot be reliably distinguished so we have to accept some loss in
# sensitivity to get a reasonably low rate of spurious OTUs.

# Analyzing a control sample such as a staggered mock community can be used to calibrate parameters to get the highest sensitivity with
# an acceptable rate of spurious OTUs.

# This page describes how you can investigate coverage indepedently of control samples, e.g. because you didn't sequence any.

# Classifying reads
# Reads can be divided into five groups as shown in the figure.

# (a) Reads that map to OTUs. If we assume that the OTUs are correct biological sequences, then these reads are correct or within 3% of
# a correct sequence. These are the reads which are counted in the OTU table if you follow my recommendations for mapping reads to OTUs.
# I think it is reasonable to consider these reads explained by the OTUs, while acknowledging that there may be some cases where multiple
# species are lumped into a single OTU.

# (b) Reads of an OTU with >3% errors. The vast majority of these will have been discarded because they do not map to an OTU.

# (c) Correct reads of chimeric amplicons. Most likely, only a small minority of OTU sequences will be chimeric because the cluster_otus
# and unoise3 commands have very effective chimera filters. If a chimera is <3% diverged from an OTU, some or all of its reads will be
# mapped to that OTU. I consider this to be harmless for all practical purposes because chimeras are much less abundant than their
# parents, so the number of chimeric reads added to the OTU count will be small compared to other fluctuations. The large majority, or
# perhaps all, chimeras in your reads will have parents that are in the OTUs because chimeras are strongly biased to form between
# high-abundance parents, for pretty obvious reasons.

# (d) Chimeric reads with errors. If a chimeric read has errors, say >3% compared with the true chimeric amplicon sequence, then it is
# impossible to distinguish it reliably from a noisy read of one of the parent sequences or a biological sequence missing from the OTUs.

# (e) Reads of biological sequences that are not in the OTUs. We would like to estimate the size of this set, but this cannot be done
# reliably because the only way to do that is to estimate the sizes of (a) through (d) and see what is left over, but the uncertainties
# in those sizes are comparable with the size of (e). The calculation ends up being something like (e) = 100% minus (90% +/- 10%),
# so you really don't have much of an idea how big (e) is -- it could be anywhere between 0% and 20%.

# Accounting for reads which do not map to an OTU at 97% identity
# Make a FASTQ file for the unmapped reads by using the -notmatchedfq option of the otutab command.

# We can safely assume that a large majority of low-quality reads which did not map belong to the OTUs and chimeras formed from the OTUs.
# This is because the OTUs represent the most abundant template sequences, and most bad reads will be derived from those templates just
# as most good reads are. With this in mind, pick an expected error threshold, say 2, and consider the reads above the threshold to be
# explained as noisy reads of known OTUs.

# usearch -fastq_filter unmapped.fq -fastq_maxee 2.0 -fastaout unmapped_hiqual.fa \
#   -fastaout_discarded unmapped_loqual.fa

# Now let's try to find chimeric reads. Use the -chimeras option of unoise3 to get the predicted chimeric amplicons. Do this even if
# you're only using 97% OTUs in your analysis, because this is the best way to get an accurate set of chimera predictions.

#  If an unmapped read is within, say, 5% of a chimera or an OTU, then it is probably a noisy read with more errors than expected from
# its Phred scores. Combine the OTUs and chimeras into one database file and use this command to check.

# usearch -usearch_global unmapped_hiqual.fa -db otus_chimeras.fa -strand plus -id 0.95 \
#   -matched unmatched_noisy.fa -notmatched unmatched_hiqual_other.fa

# Now you have the number of mapped reads and counts of low quality reads and noisy reads with apparently high quality. If you add
# these up, it should account for most of your data.

# One way (the only way?) to check for valid biological sequences that are not accounted for in the OTUs is to search the remainng
# reads against a large database like SILVA, e.g.

# usearch -usearch_global unmatched_hiqual_other.fa -db silva.udb -strand both \
#   -id 0.99 -alnout unmapped_silva.aln

# I specified a high identity threshold because if you use a lower threshold, say 0.97, you may get spurious alignments from noisy
# reads. If the identity is much below 100%, it is difficult to decide..

#########################################################################
##### NOTE: HAVEN'T BEEN ABLE TO GET THE -CHIMERAS OPTION TO WORK########
#########################################################################

echo ""
figlet -f bubble Step 11
figlet -f small Quality controlling OTUs and ZOTUs
echo ""

echo ""
figlet -f digital CREATING QUALITY CONTROL DIRECTORIES
echo ""

mkdir $quality_control_otus
mkdir $quality_control_zotus

echo ""
figlet -f digital QUALITY FILTERING UNMAPPED READS
echo ""

$usearch10 -fastq_filter $otu_table/9_unmapped.fq -fastq_maxee 2.0 -fastaout $quality_control_otus/11a_otu_unmapped_hiqual.fa \
-fastaout_discarded $quality_control_otus/11b_otu_unmapped_loqual.fa

$usearch10 -fastq_filter $zotu_table/10_unmapped.fq -fastq_maxee 2.0 -fastaout $quality_control_zotus/11a_zotu_unmapped_hiqual.fa \
-fastaout_discarded $quality_control_zotus/11b_zotu_unmapped_loqual.fa

echo ""
figlet -f digital MAPPING HIGH QUALITY UNMAPPED READS TO SILVA DATABASE
echo ""

$usearch9 -usearch_global $quality_control_otus/11a_otu_unmapped_hiqual.fa -db ZZ_Silva.SSU.udb -strand both -id 0.99 \
-alnout $quality_control_otus/11c_unmapped_silva_otus.aln

$usearch9 -usearch_global $quality_control_zotus/11a_zotu_unmapped_hiqual.fa -db ZZ_Silva.SSU.udb -strand both -id 0.99 \
-alnout $quality_control_zotus/11c_unmapped_silva_zotus.aln

echo ""
figlet -f digital OUTPUT FILES:
echo $quality_control_otus
echo $quality_control_zotus
echo 11a_otu_unmapped_hiqual.fa
echo 11b_otu_unmapped_loqual.fa
echo 11a_zotu_unmapped_hiqual.fa
echo 11b_zotu_unmapped_loqual.fa
echo 11c_unmapped_silva_otus.aln
echo 11c_unmapped_silva_zotus.aln
echo ""

# ------------------------------------------------------CORRECTING CROSS-TALK-------------------------------------------------------

# uncross command
# Commands > OTU analysis
# See also
#   Cross-talk
#   UNCROSS algorithm
#   UNCROSS paper
#   OTU table
#   Example OTU tables and reports

# The uncross command detects and filters cross-talk (sample mis-assignment) in a OTU table using the UNCROSS algorithm. In a typical
# run, about 2% of reads are assigned to the wrong sample. If some samples contain large numbers of reads for a given OTU, these often
# "bleed" into other samples which may not in fact contain that OTU. This can cause may spurious counts wihch should be zero, giving
# inflated estimates of richness, alpha diversity and beta diversity.

# You can clearly see cross-talk in this GAIIx example and this MiSeq example. You can use this example data to try the uncross command.

# Please note that I do not consider the UNCROSS algorithm to be a robust solution for cross-talk. The mechanism(s) causing cross-talk
# are not well understood. Many different indexing schemes are used. Cross-talk rates in your data may be quite different from the
# datasets on which UNCROSS was designed and tested, in which case the accuracy of UNCROSS on your data may be lower. Also,
# cross-talk may be hard or impossible to detect when the number of multiplexed samples is large, say around 100 or more. It is much
# better to use multiplexing strategies that are designed to reduce cross-talk. UNCROSS is best understood as a simplisitc hack that is
# the best we can do with exisitng data.

# Input is an OTU table in QIIME classic format generated from the all of the reads in a single run. Runs should NOT be combined for
# this analysis. It is important to include ALL samples that were sequenced in the same run, even if they contain samples for different
# experiments.

# If the run has mock community samples, mock sample names should start with "mock" (case-insensitive), e.g. Mock1, mock or mock_13.
# OTU identifiers for sequences that are in the designed mock community should contain one of the following strings (case-sensitive):
# ";mock=yes;", ";annot=perfect;" or ";annot=noisy;" The annot command can be used to generate these annotations in the sequence labels
# before generating the OTU table.

# The -tabbedout option specifies an output file in tabbed text format.

# The -report option specifies a text file name for a summary report.

# The -otutabout option specifies a filename to store the filtered OTU table. By default, entries predicted to be spurious due to
# cross-talk are set to zero and undetermined entries are kept. Specify -uncross_undet_zero to set undetermined entries to zero.

# The following options specify user-settable parameters:

# -uncross_maxxt (default 0.5). Maximum cross-talk frequency as a percentage.

# -uncross_minvalid (default 2.0). Minimum valid frequency as a percentage.

# -uncross_minvalidtotal (default 75.0). Minimum fraction of valid reads in an OTU as a percentage.

echo ""
figlet -f digital CORRECTING CROSS-TALK
echo ""

echo ""
figlet -f digital CREATING FILTERED OTU AND ZOTU DIRECTORIES
echo ""

filtered_otus="11c_filtered_OTUs"
filtered_zotus="11d_filtered_ZOTUs"

mkdir $filtered_otus
mkdir $filtered_zotus

$usearch10 -uncross $otu_table/9_otutab.txt -tabbedout $filtered_otus/11a_out.txt -report $filtered_otus/11b_rep.txt \
-otutabout $filtered_otus/11c_filtered_otu_table.txt

$usearch10 -uncross $zotu_table/10_zotutab.txt -tabbedout $filtered_zotus/11d_out.txt -report $filtered_zotus/11e_rep.txt \
-otutabout $filtered_zotus/11f_filtered_otu_table.txt

echo ""
figlet -f digital OUTPUT FILES:
echo $filtered_otus
echo $filtered_zotus
echo 11a_out.txt
echo 11b_rep.txt
echo 11c_filtered_otu_table.txt
echo 11d_out.txt
echo 11e_rep.txt
echo 11f_filtered_otu_table.txt
echo ""

# ----------------------------------------------------------ASSIGNING TAXONOMY----------------------------------------------------------

# sintax command
# Commands > Taxonomy

#  See also
#   SINTAX reference data downloads
#   Which taxonomy database should I use?
#   SINTAX algorithm
#   makeudb_sintax command
#   Can SINTAX predict species?
#   sintax_summary command

# The sintax command uses the SINTAX algorithm to predict taxonomy for query sequences
#  in FASTA or FASTQ format. See SINTAX paper for details.

# You can use the sintax_summary command to get a tabbed text file for making figures.
# The search database must have taxonomy annotations. The makeudb_sintax command can be used to create a UDB
# database, which is faster to load. See SINTAX downloads page for available reference files in FASTA format.
# See also Which database should I use?

# Taxonomy predictions are written to the -tabbedout file. The first three fields are (1) query sequence label, (2)
# prediction with boostrap values and (3) strand. If the -sintax_cutoff option is given then predictions are written a second
# time after applying the confidence threshold, keeping only ranks with high enough confidence. On V4 reads, using a cutoff of
# 0.8 gives predictions with similar accuracy to RDP at 80% bootstrap cutoff.

# The strand option must be specified.

# Multithreading is supported.

# Example

# usearch -sintax reads.fastq -db 16s.udb -tabbedout reads.sintax -strand both -sintax_cutoff 0.8

echo ""
figlet -f bubble Step 12
figlet -f small Assigning taxonomy
echo ""

echo ""
figlet -f digital CREATING TAXONOMY DIRECTORIES
echo ""

mkdir $assigned_taxonomy_otus
mkdir $assigned_taxonomy_zotus

echo ""
figlet -f digital ASSIGNING TAXONOMY TO OTUs AND ZOTUs
echo ""

$usearch9 -sintax $OTUs/otus.fa -db ZZ_Silva.SSU.udb -tabbedout $assigned_taxonomy_otus/12a_assigned_taxonomy_otus.sintax -strand both \
-sintax_cutoff 0.8

$usearch9 -sintax $ZOTUs/zotus.fa -db ZZ_Silva.SSU.udb -tabbedout $assigned_taxonomy_zotus/12b_assigned_taxonomy_zotus.sintax \
-strand both -sintax_cutoff 0.8

echo ""
figlet -f digital OUTPUT FILES:
echo $assigned_taxonomy_otus
echo $assigned_taxonomy_zotus
echo 12a_assigned_taxonomy_otus.sintax
echo 12b_assigned_taxonomy_zotus.sintax
echo ""

# -------------------------------------------------------MEASURING ALPHA DIVERSITY-------------------------------------------------------

# alpha_div command
# Commands > OTU analysis, Diversity analysis
# See also
#   alpha diversity metrics
#   beta_div command

# The alpha_div command calculates one or more alpha diversity metrics for each sample in an OTU table. The OTU table must be in QIIME
# classic format.

# The -metrics option specifies one or more metric names separated by commas. Default is all metrics. See alpha diversity metrics for
# supported metrics.

# Output is written to a table specified by the -output option. The file is in tabbed text format. Rows are samples, columns are metrics.
# This file can be loaded into spreadsheet software such as Excel for further analysis.

# Example: calculate all metrics

# usearch -alpha_div otutable.txt -output alpha.txt

# Example: calculate Gini-Simpson index

# usearch -alpha_div otutable.txt -output gini.txt -metrics gini_simpson

# Example: calculate Chao1 and Berger-Parker indexes

# usearch -alpha_div otutable.txt -output alpha.txt -metrics chao1,berger_parker

# alpha_div_rare command
# Commands > OTU analysis, Diversity analysis
# See also
#   Rarefaction
#   Abundance rarefaction

# The alpha_div_rare command calculates a rarefaction curve for an alpha diversity metric.
# Note that this command currently does not currently account for discarding singletons. I doubt it matters in practice, because other
# sources of error are probably more important, so rarefaction analysis has dubious value for marker gene OTUs. If you are interested
# in a rarefaction command which accounts for singletons, please let me know and I will implement the feature.

# Input is an OTU table in QIIME classic format.

# The -metric option specifies a metric name. Default is richness. See alpha diversity metrics for supported metrics.

# The specified metric is calculated by subsampling the OTU table at 1%, 2% ... 99%, 100% of the reads. The -method option specifies
# which subsampling method to use. Default is fast subsampling. See otutab_subsample command for discussion of subsampling and supported
# methods.

# The output file is specified by the -output option. The output file is a tabbed text file with a header followed by 100 lines, one per
# subsample size. The first field is the subsample size as a percentage.

# The output file can easily be loaded into a spreadsheet or other charting software for generating figures.

# Example

# usearch -alpha_div_rare otutab.txt -output rare.txt

echo ""
figlet -f bubble Step 13
figlet -f small Measuring alpha diversity
echo ""

echo ""
figlet -f digital CREATING ALPHA DIVERSITY DIRECTORY
echo ""

mkdir $alpha_diversity

echo ""
figlet -f digital CALCULATING ALPHA RAREFACTION CURVES
echo ""

$usearch10 -alpha_div $otu_table/9_otutab.txt -output $alpha_diversity/13a_otu_alpha.txt

$usearch10 -alpha_div $zotu_table/10_zotutab.txt -output $alpha_diversity/13b_zotu_alpha.txt

echo ""
figlet -f digital CALCULATING ALPHA DIVERSITY
echo ""

$usearch10 -alpha_div_rare $otu_table/9_otutab.txt -output $alpha_diversity/13c_otu_rare.txt

$usearch10 -alpha_div_rare $zotu_table/10_zotutab.txt -output $alpha_diversity/13d_zotu_rare.txt

echo ""
figlet -f digital OUTPUT FILES:
echo $alpha_diversity
echo 13a_otu_alpha.txt
echo 13b_zotu_alpha.txt
echo 13c_otu_rare.txt
echo 13d_zotu_rare.txt
echo ""

# -------------------------------------------------------MEASURING BETA DIVERSITY-------------------------------------------------------

# beta_div command
# Commands > OTU analysis, Diversity analysis
# See also
#   Diversity analysis
#   Beta diversity
#   beta diversity metrics
#   alpha_div command

# The beta_div command calculates one or more beta diversity metrics from an OTU table. The OTU table must be in QIIME classic format.

# The -metrics option specifies one or more metric names separated by commas. Default is all supported metric names. See beta diversity
# metrics for supported names.

# For each metric, up to three output files are generated: a distance matrix in square format, a sorted distance matrix and a tree.

# Trees are written in Newick format.

# Sorted matrices and trees are generated only for metrics which support clustering; see beta diversity metrics. In a sorted distance
# matrix, samples are sorted to bring similar samples together.

# Output file names are based on the metric name. The -filename_prefix option species a string to be added at the beginning of all
# output filenames. Typically this is a directory name, which must end with a slash or backslash. The directory must exist; usearch
# will not create a new directory. Default is no prefix.

# The ‑mx_suffix, ‑sorted_mx_suffix and -tree_suffix options specify strings to be added at the end of the output filenames. Default
# values are .txt, .sorted.txt and .tree, respectively.

# The Unifrac metric requires a tree for the OTUs, which is specified using the -tree option. The file must be in Newick format. The
# tree can be generated using the -cluster_agg command using the OTU FASTA file as input. If no tree is specified, the Unifrac metric
# will not be calcualated.

# Example: calculate all supported beta metrics and write results in current directory

# usearch -beta_div otutable.txt

# Example: calculate all supported beta metrics and write output in directory /results/beta

# usearch -beta_div otutable.txt -filename_prefix /results/beta/

# Example: calculate Jaccard only and write results in current directory

# usearch -beta_div otutable.txt -metrics jaccard

# Example: calculate Jaccard and Bray-Curtis and write results in current directory

# usearch -beta_div otutable.txt -metrics jaccard,bray_curtis

echo ""
figlet -f bubble Step 14
figlet -f small Measuring beta diversity
echo ""

echo ""
figlet -f digital CREATING BETA DIVERSITY DIRECTORY
echo ""

mkdir $beta_diversity_otus
mkdir $beta_diversity_zotus

echo ""
figlet -f digital CALCULATING BETA DIVERSITY
echo ""

$usearch10 -beta_div $otu_table/9_otutab.txt -filename_prefix $beta_diversity_otus/

$usearch10 -beta_div $zotu_table/10_zotutab.txt -filename_prefix $beta_diversity_zotus/

echo ""
figlet -f digital OUTPUT FILES:
echo $beta_diversity_otus
echo $beta_diversity_zotus
echo ""

# -------------------------------------------------BUILDING TREES AND DISTANCE MATRICES-------------------------------------------------

# cluster_agg command
# Commands > Clustering
# See also
#   cluster_otus
#   cluster_fast
#   cluster_smallmem
#   cluster_aggd

# The cluster_agg command performs agglomerative clustering of sequences in a FASTA or FASTQ file.

# Cluster linkage is specified using the ‑linkage option, which may be set to max (the default), min or avg.

# Output is reported as a tree in Newick format specified by the -treeout option and/or as a clusters file specified by the -clusterout
# option. If a clusters file is specified, then the -id option must be given to specify the identity threshold.

# The first step in the algorithm is to compute a distance matrix, which can be saved to a file by specifying the -distmxout option.
# See calc_distmx for options related to the distance matrix calculation and output file.

# calc_distmx command
# Commands > Trees and distance matrixes
# Generate a distance matrix from an input file in FASTA or FASTQ format. See also calc_distmx_smallmem.

# The distance matrix filename is specified by the -distmxout option.

# The matrix format is specified by the -format option, which can be tabbed_pairs (default), square, phylip_square or
# phylip_lower_triangle. See distance matrix for details for these file formats.

# Distance values are in the range zero (identical sequences) to one (no similarity) corresponding to the range 100% identity to 0%
# identity.

# Multithreading is supported.

# Clusters can be generated from a distance matrix with the cluster_aggd command.

# By default, pairs are prioritized by the U-sort heuristic as used in the USEARCH algorithm. This means that pairs are considered in
# decreasing order of the number of unique words (U) they have in common. Since U correlates with identity, this means that pairs are
# considered in approximately increasing order of distance. U-sorting can be turned off using the -nousort option. U-sorting plus
# additional heuristics used to find HSPs can all be disabled using the -distmx_brute option, which forces all pairs of sequences to
# be aligned. This is guaranteed to give a complete matrix, but can be much slower for large datasets. Note that low-identity pairs
# generally have little effect on clustering or tree topology, so the additional "accuracy" of a brute force calculation often has
# little biological value.

# The -sparsemx_minid option gives the minimum identity which should be written to a matrix in tabbed_pairs format.

# An identity threshold for terminating the calculation can be specified using the termid option, which is in the range 0.0 to 1.0,
# where 1.0 means identical sequences (100% sequence id). This is a speed optimization that saves time by skipping alignments of
# low-identity pairs. If a pair is encountered with fractional identity < termid, the calculation is stopped. Because U-sorted order
# does not correlate perfectly with identity, you should set termid somewhat lower than the minimum identity that you care about. For
# example, if you want all pairs with >80% id to appear in the matrix, then you might set -termid 0.7. Tests on small datasets can be
# used to tune -termid to a reasonable value. By default, termid is set to 0 and the calculation continues for all pairs that have at
# least one word in common. The word length is set by the wordlength option.

# Examples

# usearch -calc_distmx seqs.fa -distmxout mx.txt -sparsemx_minid 0.8 -termid 0.7

# usearch -calc_distmx seqs.fa -distmxout dist.tree -format phylip_lower_triangular


# cluster_aggd command
# Commands > Clustering
# The cluster_aggd command performs agglomerative clustering of a distance matrix in tabbed pairs format. To cluster sequences, use
# cluster_agg.

# The calc_distmx command can be used to calculate a distance matrix for sequences.

# The beta_div command can be used to calcualte a distance matrix for samples in an OTU table.

# Cluster linkage is specified using the ‑linkage option, which may be set to max (the default), min or avg.

# Output is reported as a tree in Newick format specified by the -treeout option and/or as a clusters file specified by the -clusterout
# option. If a clusters file is specified, then the -id option must be given to specify the identity threshold.

# Example

# usearch -cluster_aggd mx.txt -treeout clusters.tree -clusterout clusters.txt -id 0.80 -linkage min

##################################################################################################################################################
# CAN'T GET -clusterout $trees_otus/15b_otu_clusters.txt -id 0.80 -linkage min TO WORK, GETTING FOLLOWING ERROR: "TreeToClusters not implemented"#
##################################################################################################################################################

echo ""
figlet -f bubble Step 15
figlet -f small Building trees and distance matrices
echo ""

echo ""
figlet -f digital CREATING DIRECTORIES FOR TREES AND DISTANCE MATRICES
echo ""

mkdir $trees_otus
mkdir $trees_zotus

echo ""
figlet -f digital BUILDING TREES FOR OTUs AND ZOTUs
echo ""

$usearch10 -cluster_agg $OTUs/otus.fa -treeout $trees_otus/15a_otu_seqs.tree -id 0.80 -linkage min

$usearch10 -cluster_agg $ZOTUs/zotus.fa -treeout $trees_zotus/15a_zotu_seqs.tree -id 0.80 -linkage min

echo ""
figlet -f digital BUILDING DISTANCE MATRICES
echo ""

$usearch10 -calc_distmx $OTUs/otus.fa -distmxout $trees_otus/15b_otu_distance_matrix.txt -sparsemx_minid 0.8 -termid 0.7 \
-format tabbed_pairs

$usearch10 -calc_distmx $ZOTUs/zotus.fa -distmxout $trees_zotus/15b_zotu_distance_matrix.txt -sparsemx_minid 0.8 -termid 0.7 \
-format tabbed_pairs

echo ""
figlet -f digital OUTPUT FILES:
echo $trees_otus
echo $trees_zotus
echo 15a_otu_seqs.tree
echo 15a_zotu_seqs.tree
echo 15b_otu_distance_matrix.txt
echo 15b_zotu_distance_matrix.txt
echo ""

echo ""
figlet 16S Analysis Complete
echo ""
