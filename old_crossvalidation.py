#import Pybrain libraries
from pybrain.datasets            import ClassificationDataSet, SupervisedDataSet
from pybrain.utilities           import percentError
from pybrain.tools.shortcuts     import buildNetwork
from pybrain.tools.validation    import CrossValidator
from pybrain.supervised.trainers import BackpropTrainer
from pybrain.structure.modules   import SoftmaxLayer
from pybrain.tools.validation    import ModuleValidator

#Import some extra needed stuff
import numpy as np
import os, sys
import pandas as pd

#Import own libraries
from biom_parser import biom_parser as bp
from collection import collection as col
from plot_validation import validation_plot
from plot_testsets import testset_plot 
import os, sys
from cross_validation_methods import *

#import argparse and pickle
import argparse
import pickle 

parser = argparse.ArgumentParser(description=\
  'This program classifies taxonomic classifications with a PyBrian neural network')
parser.add_argument('-e',type=int, required=True,\
  help="Number of epochs to train")
parser.add_argument('-n',type=int, required=True,\
  help="Number of hidden nodes")
parser.add_argument('-pv',action="store_true", default=False, required=False,\
  help="Plot training and validation error")
parser.add_argument('-pt',action="store_true", default=False, required=False,\
  help="Plot heatmap of testsets")
parser.add_argument('-s',action="store_true",required=False,\
  help="Use pickle to save network to file for future use")
parser.add_argument('-l',type=str,required=False,\
  help="Load a trained network saved with pickle")
parser.add_argument('--train',action="store_true",required=False,\
  help="Input firectory for trainingsets")
parser.add_argument('--test',action="store_true",required=False,\
  help="Input directory for test sets")
args = parser.parse_args()


traindir = "train"
testdir = "test"
upper_taxa_level = 2
upper_taxa_level = 8
lower_taxa_level = 9
# The Columns array decides which categories to use as output nodes in the ANN. May be multiple columns. Indexed by 0.
columns = ['10']
#columns = ['9']
collection = col()

#Parse training-data
for file in os.listdir(traindir):
	if file.endswith('3'):
		a = bp(os.path.join(traindir, file), os.path.join(traindir, 'table.tsv'))
		a.createMgrast()
		a.parseMgrastCategories(columns, '2')

		#Level 1 = Root, Level 2 = Kingdom, Level 3 = Phylum...etc
		#a.set_lowest_taxa_level(lower_taxa_level)
		#a.delete_upper_taxa_level(upper_taxa_level)
		#a.convertToRelative()
		#a.convertToShares()
		a.convertToPresence()
		for otu in a.getOtuList():
			collection.addTaxa(otu)
		for categorylist in a.getCategoryList():
			for category in categorylist:
				collection.addCategory(category)
		collection.addFile(a, 'train')
#
#Parse test-data
#for file in os.listdir(testdir):
#    if file.endswith('fasta'):
#		a = bp(os.path.join(testdir,file), os.path.join(testdir, 'metadata.csv'))
#		a.createMgrast()
#		a.parseMgrastCategories('3', '1')
	        #Level 1 = Root, Level 2 = Kingdom, Level 3 = Phylum...etc
		#a.set_lowest_taxa_level(lower_taxa_level)
	        #a.delete_upper_taxa_level(upper_taxa_level)
	        #a.createDataDict()
		#a.convertToRelative()
		#a.convertToShares()
#		a.convertToPresence()
#		for otu in a.getOtuList():
#			collection.addTaxa(otu)
#		for categorylist in a.getCategoryList():
#			for category in categorylist:
#				collection.addCategory(category)
#		collection.addFile(a, 'test')

# Calculate and print number of total inputnodes (unique taxa) and total output nodes (uniqe categories to classify)
collection.setUniqueTaxa()
collection.setUniqueCategories()
print 'Unique Taxa (#input nodes): '+str(len(collection.getUniqueTaxa()))
print 'Unique Categories (#output nodes): '+str(len(collection.getUniqueCategories()))

# Create trainingsets and test sets
trainingset = collection.createAnnTrainingsets()
#print trainingset	

# Map trainingsets and test sets to PyBrain
DS = ClassificationDataSet(trainingset['input_dimension'], trainingset['output_dimension'])
for i in range(0, len(trainingset['input_arrays'])):
	DS.appendLinked(trainingset['input_arrays'][i] , trainingset['output_arrays'][i])


#Train network and plot
#if not args.l:
#        errorlist = trainer.trainUntilConvergence(maxEpochs=args.e,validationProportion=0.2)
#        if args.pv:
#                validation_plot(errorlist)

# Split dataset into five parts (1/5 test and 4/5 trainingsets)
num_splits = 5
dividend = (DS.getLength()/num_splits)
indicies = np.random.permutation(DS.getLength())
sliced_indicies = slice_list(indicies, num_splits)
trainingsets = []
testsets = []
target_actual = []
target_pred = []

#Create 5 parts, append to trainingsets
for part in sliced_indicies:
	testset = SupervisedDataSet(inp=DS['input'][part].copy(), target=DS['target'][part].copy())
	testsets.append(testset)
	trainingset = set(indicies) - set(part)
	trainingset = SupervisedDataSet(inp=DS['input'][list(trainingset)].copy(), target=DS['target'][list(trainingset)].copy())
	trainingsets.append(trainingset)


#Create ANN, run training and validation on each part
for i in range(num_splits):
	print "Part"+str(i+1)
	fnn = buildNetwork( DS.indim, args.n, DS.outdim, outclass=SoftmaxLayer, fast=False )
	trainer = BackpropTrainer(fnn, dataset=trainingsets[i], momentum=0.01, verbose=True, weightdecay=0.0001)
	trainer.trainEpochs(args.e)
	out = fnn.activateOnDataset(testsets[i])
#	print list(collection.unique_categories)
#	print testsets[i]
	actual = testsets[i]['target'].argmax(axis=1)
#	print actual.tolist()
	pred = out.argmax(axis=1)
#	print pred.tolist()
	new_actual = [list(collection.unique_categories)[x] for x in actual.tolist()]
	new_pred = [list(collection.unique_categories)[x] for x in pred.tolist()]
#	print new_actual
#	print new_pred
	target_actual += new_actual
	target_pred += new_pred


y_actu = pd.Series(target_actual ,name='actual')
y_pred = pd.Series(target_pred ,name='predicted')
cm_margins = pd.crosstab(y_actu, y_pred, rownames=['Actual'], colnames=['Predicted'], margins=True)
cm = pd.crosstab(y_actu, y_pred, rownames=['Actual'], colnames=['Predicted'])
cm_norm = pd.crosstab(y_actu, y_pred, rownames=['Actual'], colnames=['Predicted'])
cm_norm = cm_norm / cm_norm.sum(axis=1)
print cm_margins
print "Number of target sets: "+str(DS.getLength())

True_positives = [int('1') for x in range(len(y_actu.tolist())) if y_actu.tolist()[x] == y_pred.tolist()[x]]
print "Accuracy: "+str(float(sum(True_positives)) / float(cm_margins['All'][-1]))
#print "Accuracy: "+str(float(np.diag(cm).sum()) / float(cm_margins['All'][-1]))

exit()


#Load network
if args.l:
        f = open (args.l, "r")
        fnn = pickle.load(f)

#Save network 
if args.s:
        with open("network.obj", "w") as f:
                pickle.dump(fnn, f)
# Plot heatmap
if args.pt:
        xcategories = list(collection.unique_categories)
        ytestsets = []
        zdata = []

# Apply to testsets
testsets = {}
for testset in collection.test.keys():
        anntestset = collection.createAnnTestset(testset)
	DStest = ClassificationDataSet(trainingset['input_dimension'], trainingset['output_dimension'])
	DStest.appendLinked(anntestset['input_arrays'][0], anntestset['output_arrays'][0])
#	for i in range(0, len(anntestset['input_arrays'])):
#		DStest.appendLinked(anntestset['input_arrays'][i], anntestset['output_arrays'][i])
	print testset
	out = fnn.activateOnDataset(DStest)
	print out[0].tolist()
	print list(collection.unique_categories)
	sortedout = sorted(zip (out[0].tolist(),list(collection.unique_categories)))
	#Plot Heatmap
	if args.pt:
		tmp1 = testset.split('/')
		tmp2 = tmp1[-1].split('.')
		testset = tmp2[0]
		ytestsets.append(testset)
		zdata.append(out[0].tolist())
	#sortedout = [x for (y,x) in sortedout]
	#print sortedout
#	print collection.unique_categories
#	print out

#Plot heatmap
if args.pt:
        testset_plot(xcategories, ytestsets, zdata)

