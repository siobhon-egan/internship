#!/bin/bash

# Enter FWD primer sequence 5'-3' (degenerate bases OK)
fwd_primer="AGAGTTTGATCCTGGCTYAG"
# Enter REV primer sequence 5'-3' (degenerate bases OK)
rev_primer="TGCTGCCTCCCGTAGGAGT"

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Triming primers and distal bases
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo ""

# At the moment this usearch command can only take .fasta as input so can only be done
# after QF

# Creating working directories

mkdir trimmed_data_test
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

for file4 in quality_filered_test/*.fasta
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Trimming primers step 1: finding seqs with FWD primer
        echo input is:
        echo ${file4}

    usearch9.2 -search_oligodb ${file4} -db fwd_primer_db.fasta -strand both -matched "seqs_w_fwd_primer/$(basename ${file4})" -notmatched "seqs_wo_fwd_primer/$(basename ${file4})"
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

    usearch9.2 -search_oligodb ${file5} -db rev_primer_db.fasta -strand both -matched "seqs_w_fwd_and_rev_primer/$(basename ${file5})" -notmatched "seqs_w_fwd_butnot_rev_primer/$(basename ${file5})"
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

    usearch8 -search_pcr ${file6} -db both_primers_db.fasta -strand both -maxdiffs ${pcr_missmatches} -pcr_strip_primers -ampout "trimmed_data_test/$(basename ${file6} .fasta).fasta"
done
