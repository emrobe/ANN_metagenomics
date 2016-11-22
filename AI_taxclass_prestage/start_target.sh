#!/bin/sh

#$ -N AI_taxclass_prestage
#$ -S /bin/bash
#$ -o /home/emr023/AI_taxclass_prestage/start.log
#$ -j y

# Table.tsv needs an added first column with ID's starting from one (1,2,3,4,Xn...) which will be referenced by ${SGE_TASK_ID}
cd /home/emr023/AI_taxclass_prestage
#/share/apps/python/Python-2.7.8/python 
#source 
/home/emr023/virtualenv-13.1.2/python2.7/bin/python2.7 create_ai_targetsets_mgrast.py -n 8 ${SGE_TASK_ID} 

