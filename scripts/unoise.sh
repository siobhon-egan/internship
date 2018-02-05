echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo UNOISE on all
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir ${unoise_all}

cat ${SF}/*.fasta > ${unoise_all}/all_samples_SF.fasta

cd ${SF}
gzip *.fasta
cd ..

cd ${SF_derep}
gzip *.fasta
cd ..

cd ${unoise_all}

usearch9.2_linux -fastx_uniques all_samples_SF.fasta -fastaout all_samples_SF_DR.fasta -sizeout

gzip all_samples_SF.fasta
