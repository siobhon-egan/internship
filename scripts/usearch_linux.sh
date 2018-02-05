#!/bin/bash
#
#	Requirements: linux Nimbus VM requires usearch9.2 and usearch10 (ensure they are linux downloads)
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

# Usearch 8 and 9 files as it is on the $PATH
usearch9="usearch9.2_linux"
usearch8="usearch8.0_linux32"
# Enter raw data directorry
raw_data="/data/subset/raw_data"
# Enter directory for merged output
read_summary="0.read_summary"
# Enter directory for merged output
merged_data="1.merged_sequences"
# Enter minimum merge overlap - 50 bp minimum
overlap="50"
# Enter directory for quality filtered output
QF="quality_filtered"
# Enter % expected error rate threshold for seq quality filtering (1% = 0.01)
error_rate="0.01"
# Enter min length of sequence for trimming in bp (eg. to keep all seqs above 200 bp enter "200")
minlen="150"
# Enter directory for trimmed data
trimmed_data="trimmed_data"
# Enter FWD primer sequence 5'-3' (degenerate bases OK)
fwd_primer="AGAGTTTGATCCTGGCTYAG"
# Enter REV primer sequence 5'-3' (degenerate bases OK)
rev_primer="TGCTGCCTCCCGTAGGAGT"
# Enter number of primer missmatches allowed
pcr_missmatches="0"
# Enter directory for labeled data
labeled_data="2.labeled_data"
# Enter directory for derep seq (unique seqs)
derep_dir="3.dereplicated_sequences"
# Enter directory for singleton seqs
low_abund_seqs="4.singleton_sequences"
# Enter directory for singleton filtered data
SF="5.singleton_filtered"
# Enter max replicate cluster size (eg. to remove singletons enter 1, for duplicates enter 2)
maxsize="1"
# Enter directory for SF-derep data
SF_derep="6.dereplicated_filtered_data"

##########################################################################################
# DO NOT EDIT BELOW THIS LINE
##########################################################################################
mkdir $read_summary

 for fq in $raw_data/*R1*.fastq
 do
  $usearch9 -fastx_info $fq -output $read_summary/1a_fwd_fastq_info.txt
 done

 for fq in $raw_data/*R2*.fastq
 do
  $usearch9 -fastx_info $fq -output $read_summary/1b_rev_fastq_info.txt
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

	$usearch9 -fastq_mergepairs ${file1} -reverse "${raw_data}/$(basename -s R1_001.fastq ${file1})R2_001.fastq" -fastqout "working1/$(basename "$file1")" -fastq_minovlen ${overlap} -report ${merged_data}/fastq_mergepairs_report_txt
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

# Zip fastq files in ${raw_data}
cd ${raw_data}
gzip *.fastq
cd ..

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

	$usearch9 -fastq_filter ${file3} -fastaout "${QF}/$(basename "$file3" .fastq).fasta" -fastq_maxee_rate ${error_rate} -fastq_minlen ${minlen}
done

#Zip fastq files in ${merged_data}
cd ${merged_data}
gzip *.fastq
cd ..
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
		mkdir seqs_w_fwd_primer
		mkdir seqs_wo_fwd_primer
		mkdir seqs_w_fwd_and_rev_primer
		mkdir seqs_w_fwd_butnot_rev_primer

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

	$usearch9 -search_oligodb ${file4} -db fwd_primer_db.fasta -strand both -matched "seqs_w_fwd_primer/$(basename ${file4})" -notmatched "seqs_wo_fwd_primer/$(basename ${file4})"
done
#*****************************************************************************************
# Step 2: Finding seqs with FWD and REV primers

for file5 in seqs_w_fwd_primer/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Trimming primers step 2: finding seqs with FWD and REV primer
		echo input is:
		echo ${file5}

	$usearch9 -search_oligodb ${file5} -db rev_primer_db.fasta -strand both -matched "seqs_w_fwd_and_rev_primer/$(basename ${file5})" -notmatched "seqs_w_fwd_butnot_rev_primer/$(basename ${file5})"
done
#*****************************************************************************************
# Step 3: Trimming FWD and REV primers

for file6 in seqs_w_fwd_and_rev_primer/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Trimming primers step 3: removing FWD and REV primers
		echo input is:
		echo ${file6}

	$usearch8 -search_pcr ${file6} -db both_primers_db.fasta -strand both -maxdiffs ${pcr_missmatches} -pcr_strip_primers -ampout "${trimmed_data}/$(basename ${file6} .fasta).fasta"
done

gzip *.fasta
#*****************************************************************************************
# Removing working directories

#		rm -r seqs_w_fwd_primer
#		rm -r seqs_wo_fwd_primer
#		rm -r seqs_w_fwd_and_rev_primer
#		rm -r seqs_w_fwd_butnot_rev_primer

# Zip working directories
cd seqs_w_fwd_primer
gzip *.fasta
cd ..
cd seqs_wo_fwd_primer
gzip *.fasta
cd ..
cd seqs_w_fwd_and_rev_primer
gzip *.fasta
cd ..
cd seqs_w_fwd_butnot_rev_primer
gzip *.fasta
cd ..

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
		rm -r ${QF}
		rm -r ${trimmed_data}

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

	$usearch9 -fastx_uniques ${file9} -fastaout "${derep_dir}/$(basename "$file9" .fasta)_DR.fasta" -sizeout
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

	$usearch9 -sortbysize ${file10} -fastaout "${low_abund_seqs}/$(basename "$file10" .fasta).fasta" -maxsize ${maxsize}
done

cd ${derep_dir}
gzip *.fasta
cd ..
#*****************************************************************************************
# Step 3: Mapping reads

for file11 in ${labeled_data}/*.fasta
	do

		echo ""
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo Removing singletons step 3: mapping reads to low abundant uniques
		echo input is:
		echo ${file11}

	$usearch9 -search_exact ${file11} -db "${low_abund_seqs}/$(basename "$file11" .fasta)_DR.fasta" -strand plus -notmatched "${SF}/$(basename "$file11" .fasta).fasta"
done

# Zip fasta files in ${labeled_data}

cd ${labeled_data}
gzip *.fasta
cd ..

cd ${low_abund_seqs}
gzip *.fasta
cd ..
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

	$usearch9 -fastx_uniques ${file15} -fastaout "${SF_derep}/$(basename "$file15" .fasta).fasta" -sizeout
done
