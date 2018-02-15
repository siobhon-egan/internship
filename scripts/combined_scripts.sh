#!/bin/bash
#
#	Requirements: usearch9.2 and usearch10 must be installed on the PATH as "usearch9.2" and "usearch10" - it does not matter whether
#	you have the 32-bit of 64-bit version.
#	I have a 64 bit version of usearch9.2, so I use that for steps that need higher RAM requiremnets.
#	This script will work in unix and linux environments.
#
#	This script taked raw MiSeq demultiplexed .fastq files for input and performs the following tasks:
#
#	1) Merging of illumina paired reads
#	2) Quality filtering of sequence data and removal of short dimer seqs
#	3) Trimming of primer sequences and distal bases and removal of sequences without correct primer sequences
#	4) Renaming files with USEARCH labels "barcodelabel=sample_id;sequence_id"
#	5) Removing low abundant sequences
#	6) OTU clusting and chimera detection with UNOISE3
##########################################################################################
#	Input raw unmerged filenames must be named "sample_id_SXXX_L001_R1_001.fastq" (read 1)
#	and "sample_id_SXXX_L001_R2_001.fastq" (read 2) and all deposited in a directory specified
#	by the "$raw_data" variable. "SXXX" is the sample number given by the MiSeq
#	in the order you entered them in the Sample Sheet.
#
#	Before use: $chmod 775 this_script.sh
#	To run: $./this_script.sh
#	This script will read any input directory specified by the "raw_data" variable, but will deposit all
#	output into the current working diretory.
##########################################################################################

# Enter raw data directorry
raw_data="raw_data"
# Enter raw data directorry
read_summary="1.read_sumary"
# Enter directory for merged output
merged_data="2.merged_data"
# Enter max diff for merging - default 5 but should increase for paired end
maxdiffs="15"
# Enter minimum merge overlap - 50 bp minimum
overlap="50"
# Enter directory for quality filtered output
QF="3.quality_filtered"
# Enter % expected error rate threchold for seq quality filtering (1% = 0.01)
error_rate="0.01"
# Enter min length of sequence for trimming in bp (eg. to keep all seqs above 200 bp enter "200")
minlen="150"
# Enter directories for trimmed data
trimmed_data="4.trimmed_data"
seqs_w_fwd_primer="5.seqs_w_fwd_primer"
seqs_wo_fwd_primer="6.seqs_wo_fwd_primer"
seqs_w_fwd_and_rev_primer="7.seqs_w_fwd_and_rev_primer"
seqs_w_fwd_butnot_rev_primer="8.seqs_w_fwd_butnot_rev_primer"
# Enter FWD primer sequence 5'-3' (degenerate bases OK)
fwd_primer="AGAGTTTGATCCTGGCTYAG"
# Enter REV primer sequence 5'-3' (degenerate bases OK)
rev_primer="TGCTGCCTCCCGTAGGAGT"
# Enter number of primer missmatches allowed
pcr_missmatches="0"
# Enter directory for labeled data
labeled_data="9.labelled_data"
# Enter directory for derep seq (unique seqs)
derep_dir="10.derep_unique_seqs"
# Enter directory for singleton filtered data
SF="11.singleton_filtered"
# Enter max replicate cluster size (eg. to remove singletons enter 1, for duplicates enter 2)
maxsize="1"
# Enter directory for SF-derep data
SF_derep="12.single_filtered_derep_data"
# Enter directory for singleton seqs
low_abund_seqs="13.singleton_sequences"
# Enter cluster directory
cluster="14.clustered_seqs"

##########################################################################################
# DO NOT EDIT BELOW THIS LINE
##########################################################################################

# Unzip raw fastq.gz files

cd ${raw_data}
gunzip *.fastq.gz
cd ..

##########################################################################################

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Read summary
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir ${read_summary}

 for fq in $raw_data/*R1*.fastq
 do
  usearch10 -fastx_info $fq -output $read_summary/1a_fwd_fastq_info.txt
 done

 for fq in $raw_data/*R2*.fastq
 do
  usearch10 -fastx_info $fq -output $read_summary/1b_rev_fastq_info.txt
 done

##########################################################################################

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Merging paried illumina sequences
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		mkdir ${merged_data}
		mkdir working1
#*****************************************************************************************
# Step1: merge data with usearch9 -fastq-filter

for file1 in ${raw_data}/*R1_001.fastq
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Merging paired reads
		echo forward reads are:
		echo $(basename ${file1})
		echo reverse reads are:
		echo $(basename ${file1} R1_001.fastq)R2_001.fastq

	usearch9.2 -fastq_mergepairs ${file1} -reverse "${raw_data}/$(basename -s R1_001.fastq ${file1})R2_001.fastq" -fastqout "working1/$(basename "$file1")" -fastq_maxdiffs ${maxdiffs} -fastq_minovlen ${overlap} -report ${merged_data}/2a_merging_seqs_report.txt -tabbedout ${merged_data}/2b_tabbedout.txt
done

#*****************************************************************************************
# Step 2: Remove "_L001_R1_001" from filenames

for file2 in working1/*.fastq
	do

		rename="$(basename ${file2} _L001_R1_001.fastq).fastq"

		mv ${file2} ${merged_data}/${rename}
done

#*****************************************************************************************
# Removing working directory

		rm -r working1

##########################################################################################
##########################################################################################

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Quality control and removing dimer seqs
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo ""

		mkdir ${QF}

#*****************************************************************************************
for file3 in ${merged_data}/*.fastq
	do
		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Quality control and removing dimer seqs
		echo input is:
		echo ${file3}

	usearch9.2 -fastq_filter ${file3} -fastaout "${QF}/$(basename "$file3" .fastq).fasta" -fastq_maxee_rate ${error_rate} -fastq_minlen ${minlen}
done

#*****************************************************************************************

##########################################################################################
##########################################################################################

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Triming primers and distal bases
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo ""

# At the moment this usearch command can only take .fasta as input so can only be done
# after QF

# Creating working directories

mkdir ${trimmed_data}
mkdir ${seqs_w_fwd_primer}
mkdir ${seqs_wo_fwd_primer}
mkdir ${seqs_w_fwd_and_rev_primer}
mkdir ${seqs_w_fwd_butnot_rev_primer}

# Creating FWD primer db

		echo ">fwd_primer" > fwd_primer_db.fasta
		echo ${fwd_primer} >> fwd_primer_db.fasta

# Creating REV primer db

		echo ">rev_primer" > rev_primer_db.fasta
		echo ${rev_primer} >> rev_primer_db.fasta

# Creating FWD and REV primer db

		echo ">fwd_primer" > both_primers_db.fasta
		echo ${fwd_primer} >> both_primers_db.fasta
		echo ">rev_primer" >> both_primers_db.fasta
		echo ${rev_primer} >> both_primers_db.fasta

#*****************************************************************************************
# Step 1: Finding seqs with FWD primer

for file4 in ${QF}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Trimming primers step 1: finding seqs with FWD primer
		echo input is:
		echo ${file4}

	usearch9.2 -search_oligodb ${file4} -db fwd_primer_db.fasta -strand both -matched "${seqs_w_fwd_primer}/$(basename ${file4})" -notmatched "${seqs_wo_fwd_primer}/$(basename ${file4})"
done
#*****************************************************************************************
# Step 2: Finding seqs with FWD and REV primers

for file5 in ${seqs_w_fwd_primer}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Trimming primers step 2: finding seqs with FWD and REV primer
		echo input is:
		echo ${file5}

	usearch9.2 -search_oligodb ${file5} -db rev_primer_db.fasta -strand both -matched "${seqs_w_fwd_and_rev_primer}/$(basename ${file5})" -notmatched "${seqs_w_fwd_butnot_rev_primer}/$(basename ${file5})"
done
#*****************************************************************************************
# Step 3: Trimming FWD and REV primers

for file6 in ${seqs_w_fwd_and_rev_primer}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Trimming primers step 3: removing FWD and REV primers
		echo input is:
		echo ${file6}

	usearch8 -search_pcr ${file6} -db both_primers_db.fasta -strand both -maxdiffs ${pcr_missmatches} -pcr_strip_primers -ampout "${trimmed_data}/$(basename ${file6} .fasta).fasta"
done

#*****************************************************************************************
# Removing working directories

#		rm -r seqs_w_fwd_primer
#		rm -r seqs_wo_fwd_primer
#		rm -r seqs_w_fwd_and_rev_primer
#		rm -r seqs_w_fwd_butnot_rev_primer

##########################################################################################
##########################################################################################

# For this script to run correctly input fasta label must be formatted >sequence_id and filename must be sample_id.fasta.
# Result will be ">barcodelabel=sample_id;sequenceid"


echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Renameing sequences with ">barcodelabel=sample_id;sequence_id"
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		mkdir ${labeled_data}
		mkdir working2

#*****************************************************************************************
# Step 1: Remove ">" from start of sequence_ID

for file7 in ${trimmed_data}/*.fasta
	do

		sed -e 's/>/>barcodelabel=;/g' ${file7} > working2/$(basename "$file7" .fasta).txt
done

#*****************************************************************************************
# Step 2: Add sample_ID (should be filename) to produce ">barcodelabel=sample_ID;sequence_ID"

for file8 in working2/*.txt
	do

		sample_id=$(basename ${file8} .txt)
		echo ${sample_id}

	sed -e "s/;/${sample_id};/g" ${file8} > "${labeled_data}/$(basename "$file8" .txt).fasta"
done


#*****************************************************************************************
# Remove working directories

		rm -r working2

##########################################################################################
##########################################################################################


echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Removing low abundant seqs singletons per sample
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo ""

# Creating directories

mkdir ${derep_dir}
mkdir ${SF}
mkdir ${low_abund_seqs}

#*****************************************************************************************
# Step 1: Dereplicating

for file9 in ${labeled_data}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Removing singletons step 1: derep_fulllength
		echo input is:
		echo ${file9}

	usearch9.2 -fastx_uniques ${file9} -fastaout "${derep_dir}/$(basename "$file9" .fasta).fasta" -sizeout
done


#*****************************************************************************************
# Step 2: Filtering low abundant seqs {maxsize}

for file10 in ${derep_dir}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Removing singletons step 2: sorting uniques
		echo input is:
		echo ${file10}

	usearch9.2 -sortbysize ${file10} -fastaout "${low_abund_seqs}/$(basename "$file10" .fasta).fasta" -maxsize ${maxsize}
done

#*****************************************************************************************
# Step 3: Mapping reads

for file11 in ${labeled_data}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Removing singletons step 3: mapping reads to low abundant uniques
		echo input is:
		echo ${file11}

	usearch9.2 -search_exact ${file11} -db "${low_abund_seqs}/$(basename "$file11" .fasta).fasta" -strand plus -notmatched "${SF}/$(basename "$file11" .fasta).fasta"
done


#*****************************************************************************************

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Dereplicating SF files
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		mkdir ${SF_derep}

#*****************************************************************************************
for file15 in ${SF}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Derep_fulllength SF data
		echo input is:
		echo ${file15}

	usearch9.2 -fastx_uniques ${file15} -fastaout "${SF_derep}/$(basename "$file15" .fasta).fasta" -sizeout
done


#*****************************************************************************************

##########################################################################################
##########################################################################################

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo CLUSTERING
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir ${cluster}

cat ${SF}/*.fasta > ${cluster}/all_samples_SF.fasta

cd ${cluster}

echo -----------------------------------------------------------------------------
echo OTU CLUSTERING - UPARSE
echo -----------------------------------------------------------------------------

mkdir ${uparse_otus}
cd ${uparse_otus}

usearch9.2 -fastx_uniques ../all_samples_SF.fasta -fastaout all_samples_SF_DR.fasta -sizeout

usearch9.2 -cluster_otus all_samples_SF_DR.fasta -otus uparse_otus.fasta -relabel OTU

usearch9.2 -usearch_global ../all_samples_SF.fasta -db uparse_otus.fasta -strand both -id 0.97 -otutabout uparse_otu_tab.txt -biomout uparse_otu_biom.biom

# Make OTU tree Newick format
usearch9.2 -cluster_agg uparse_otus.fasta -treeout otus.tree


cd ..


echo -----------------------------------------------------------------------------
echo ZOTU CLUSTERING - UNOISE3
echo -----------------------------------------------------------------------------

mkdir ${unoise3_zotus}
cd ${unoise3_zotus}

cp ../${uparse_otus}/all_samples_SF_DR.fasta .

usearch10 -unoise3 all_samples_SF_DR.fasta -zotus unoise_zotus.fasta -tabbedout unoise_tab.txt

usearch10 -fastx_relabel unoise_zotus.fasta -prefix Otu -fastaout unoise_zotus_relabelled.fasta -keep_annots

usearch10 -otutab all_samples_SF_DR.fasta -zotus unoise_zotus_relabelled.fasta -otutabout unoise_otu_tab.txt -biomout unoise_otu_biom.biom -mapout unoise_map.txt -notmatched unoise_notmatched.fasta -dbmatched dbmatches.fasta -sizeout

# Make ZOTU tree Newick format
usearch9.2 -cluster_agg unoise_zotus_relabelled.fasta -treeout zotus.tree
cd ..

cd ..



##########################################################################################
##########################################################################################
