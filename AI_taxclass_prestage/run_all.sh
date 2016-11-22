time python 16S-tool/filter_16s_new.py -t -n 4 -o ~/AI_taxclass_prestage/tmp.fasta a.fasta &&
time blastn -query ~/AI_taxclass_prestage/tmp.fasta -task megablast -db /home/emr023/AI_taxclass_prestage/silvamod/silvamod/silvamod.fasta -out blast_out -num_threads 4 -outfmt 5 -num_alignments 100 &&
time python LCAClassifier/classify.py blast_out -1 1 -2 2 -3 3 -p &&
rm tmp.fasta &&
rm blast_out 
