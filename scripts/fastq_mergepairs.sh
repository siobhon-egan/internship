#!/bin/bash

# This script works with usearch9.2 and usearch10
# This piece of code works on a subset of raw fastq reads to test for min merge overlap - default = 50 bp

mkdir merging_test

# Step1: merge data with usearch10 -fastq_mergepairs

for file1 in raw_data_subset/*R1_001.fastq
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Merging paired reads
        echo forward reads are:
        echo $(basename ${file1})
        echo reverse reads are:
        echo $(basename ${file1} R1_001.fastq)R2_001.fastq

    usearch10 -fastq_mergepairs ${file1} -reverse "raw_data_subset/$(basename -s R1_001.fastq ${file1})R2_001.fastq" -fastqout "merging_test/$(basename "$file1")" -fastq_minovlen 50 -report merging_test/report.txt
done

# Step 2: Remove "_L001_R1_001" from filenames

for file2 in merging_test/*.fastq
    do

        rename="$(basename ${file2} _L001_R1_001.fastq).fastq"

        mv ${file2} merging_test/${rename}
done
