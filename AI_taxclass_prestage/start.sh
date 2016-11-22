#!/bin/sh

#$ -N AI_taxclass_prestage
#$ -S /bin/bash
#$ -o AI_taxclass_prestage/start.log
#$ -j y

# Table.tsv needs an added first column with ID's starting from one (1,2,3,4,Xn...) which will be referenced by ${SGE_TASK_ID}

# cd working_dir/AI_taxclass_prestage

python create_ai_datasets_mgrast.py -n 8 ${SGE_TASK_ID} 

