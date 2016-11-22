# ANN_metagenomics
This repository contains the implementation of an artificial neural net in PyBrain for classifying metagenomic taxonomic classifications

### AI_taxclass_prestage/
Contains all scripts, bins and databases relevant for downloading and processing traingsets from MG-RAST. This implementation uses the Oracle Grid Engine.

Download and processing is wrapped by start.sh (training sets) and start_target.sh (test sets). Processing of training sets needs the reference file 'table.tsv' which references all relevant training sets with an inhouse ID (first column), which is passed to $SGE_TASK_ID for parallel computation. (qsub start.sh -t 1-X)

Processing of test sets is wrapped by start_tartget.sh. This script uses the reference file 'target_sampels/metadata.csv' which references all test sets in a similar manner. Raw sequence data for selected test sets is not included but can be downloaded via ENA accessions. These need to be in fasta format to be processed.

.results/ contains all processed training sets 

.target_samples/results contains all processed test sets

### computed/
Contains all result files relevant for cross validation and classification of selected test sets.

.crossvalidation/ contains bash files to compute, accuracy measurements and ploting scripts to produce surface plots

.30hidden200epochs/ contains heatmap, output and network.obj, a gzipped pickle object of the saved network (needs to be unpacked to use)

### cured_results/
Contains the original table from MG-RAST, a script to cure taxonomically classified samples and the resulting processed table which contain all relevant samples.

### test/ train/ train_subset/
Processed test sets and training sets respectively. 

Crossvalidation can be run with crossvalidation.py
Classification of novel samples can be run with classify.py 
