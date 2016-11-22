# ANN_metagenomics
This repository contains the implementation of an artificial neural net in PyBrain for classifying metagenomic taxonomic classifications

# computed/
Contains all result files relevant for cross validation and classification of selected test sets.
.crossvalidation/ contains bash files to compute, accuracy measurements and ploting scripts to produce surface plots
.30hidden200epochs/ contains heatmap, output and network.obj, a gzipped pickle object of the saved network (needs to be unpacked to use)

# cured_results/
Contains the original table from MG-RAST, a script to cure taxonomically classified samples and the resulting processed table which contain all relevant samples.

# test/ train/ train_subset/
Processed test sets and training sets respectively. 

Crossvalidation can be run with crossvalidation.py
Classification of novel samples can be run with classify.py 
