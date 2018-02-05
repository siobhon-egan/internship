# Usearch9.2 script to merge dataset

# Enter name of usearch9.2
usearch9.2="usearch9.2"
#Enter raw data directorry
raw_data="raw_data"
# Enter directory for merged output
merged_data="1.merged_data"
# Enter minimum merge overlap - 50 bp minimum
overlap="50"

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

	$usearch9.2 -fastq_mergepairs ${file1} -reverse "${raw_data}/$(basename -s R1_001.fastq ${file1})R2_001.fastq" -fastqout "working1/$(basename "$file1")" -fastq_minovlen ${overlap} -report ${merged_data}/1a_merging_seqs_report.txt -tabbedout ${merged_data}/1b_tabbedout.txt
done

#*****************************************************************************************
# Step 2: Remove "_L001_R1_001" from filenames

for file2 in working1/*.fastq
	do

		rename="$(basename ${file2} _L001_R1_001.fastq)_MG.fastq"

		mv ${file2} ${merged_data}/${rename}
done
#*****************************************************************************************
# Removing working directory

		rm -r working1
