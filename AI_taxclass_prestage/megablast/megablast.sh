#!/bin/sh
blastn -query tmp.fasta -task megablast -db /home/emr023/AI_taxclass_prestage/silvamod/silvamod/silvamod.fasta -out blast_out -num_threads 4 -outfmt 5 -num_alignments 100
