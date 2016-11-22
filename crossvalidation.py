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
parser.add_argument('--train',type=str,required=False,\
  help="Input firectory for trainingsets")
args = parser.parse_args()


upper_taxa_level = 2
upper_taxa_level = 8
lower_taxa_level = 9
# The Columns array decides which categories to use as output nodes in the ANN. May be multiple columns. Indexed by 0.
columns = ['10']
collection = col()

#Parse training-data
for file in os.listdir(args.train):
	if file.endswith('3'):
		a = bp(os.path.join(args.train, file), os.path.join(args.train, 'table.tsv'))
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
	actual = testsets[i]['target'].argmax(axis=1)
	pred = out.argmax(axis=1)
	new_actual = [list(collection.unique_categories)[x] for x in actual.tolist()]
	new_pred = [list(collection.unique_categories)[x] for x in pred.tolist()]
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

exit()

