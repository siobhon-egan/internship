# IDENDTIFY MISSING OTUS
# https://drive5.com/usearch/manual/otu_qc_missing.html
# Find the labels of the missing OTUs and make a FASTA file missing.fa with the sequences. You can do this by cutting the labels out of the OTU table and grepping the labels from the FASTA file with the OTU sequences, then using the Linux uniq command to find labels which appear only in the FASTA file. Finally, use the fastx_getseqs command to extract those sequences. For example,
# Prints out missing otus in a text file that you can copy and paste into spreadsheet and put all zero's against sample abundance
cut -f1 otutab.txt | grep -v "^#" > table_labels.txt
grep "^>" otus.fa | sed "-es/>//" > seq_labels.txt
sort seq_labels.txt table_labels.txt table_labels.txt | uniq -u > missing_labels.txt
usearch -fastx_getseqs otus.fa -labels missing_labels.txt -fastaout missing.fa
