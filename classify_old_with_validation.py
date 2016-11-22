from pybrain.datasets            import ClassificationDataSet
from pybrain.utilities           import percentError
from pybrain.tools.shortcuts     import buildNetwork
from pybrain.supervised.trainers import BackpropTrainer
from pybrain.structure.modules   import SoftmaxLayer

from pylab import ion, ioff, figure, draw, contourf, clf, show, hold, plot
from scipy import diag, arange, meshgrid, where
from numpy.random import multivariate_normal

from biom_parser import biom_parser as bp
from collection import collection as col
from plot_validation import validation_plot 
import os, sys
from pprint import pprint

#import arac
import argparse
import pickle 

parser = argparse.ArgumentParser(description=\
  'This program classifies taxonomic classifications with a PyBrian neural network')
parser.add_argument('-e',type=int, required=False,\
  help="Max number of epochs to train (using Train untill Convergence)")
parser.add_argument('-p',action="store_true", default=False, required=False,\
  help="Plot training and validation error")
parser.add_argument('-s',action="store_true",required=False,\
  help="Use pickle to save network to file for future use")
parser.add_argument('-l',type=str,required=False,\
  help="Load a trained network saved with pickle")
parser.add_argument('--train',action="store_true",required=False,\
  help="Input firectory for trainingsets")
parser.add_argument('--test',action="store_true",required=False,\
  help="Input directory for test sets")
args = parser.parse_args()


traindir = "pre_train"
testdir = "pre_test"
upper_taxa_level = 2
upper_taxa_level = 8
lower_taxa_level = 9
# The Columns array decides which categories to use as output nodes in the ANN. May be multiple columns. Indexed by 0.
columns = ['10']

collection = col()

#Parse training-data
for file in os.listdir(traindir):
	if file.endswith('3'):
		a = bp(os.path.join(traindir, file), os.path.join(traindir, 'table.tsv'))
		#a.createEbiBiomList()
		a.createMgrast()
		a.parseMgrastCategories(columns)

		#Level 1 = Root, Level 2 = Kingdom, Level 3 = Phylum...etc
		#a.set_lowest_taxa_level(lower_taxa_level)
		#a.delete_upper_taxa_level(upper_taxa_level)
		#a.createDataDict()
		#a.convertToRelative()
		#a.convertToShares()
		a.convertToPresence()
		for otu in a.getOtuList():
			collection.addTaxa(otu)
		for categorylist in a.getCategoryList():
			for category in categorylist:
				collection.addCategory(category)
		collection.addFile(a, 'train')

#Parse test-data
for file in os.listdir(testdir):
    if file.endswith('3'):
		a = bp(os.path.join(testdir,file), os.path.join(testdir, 'table.tsv'))
        #a.createEbiBiomList()
		a.createMgrast()
		a.parseMgrastCategories(columns)
        #Level 1 = Root, Level 2 = Kingdom, Level 3 = Phylum...etc
		#a.set_lowest_taxa_level(lower_taxa_level)
        #a.delete_upper_taxa_level(upper_taxa_level)
        #a.createDataDict()
		#a.convertToRelative()
		#a.convertToShares()
		a.convertToPresence()
		for otu in a.getOtuList():
			collection.addTaxa(otu)
		for categorylist in a.getCategoryList():
			for category in categorylist:
				collection.addCategory(category)
		collection.addFile(a, 'test')

# Calculate and print number of total inputnodes (unique taxa) and total output nodes (uniqe categories to classify)
collection.setUniqueTaxa()
collection.setUniqueCategories()
print 'Unique Taxa (#input nodes): '+str(len(collection.getUniqueTaxa()))
print 'Unique Categories (#output nodes): '+str(len(collection.getUniqueCategories()))

# Create trainingsets and test sets
trainingset = collection.createAnnTrainingsets()
testset = collection.createAnnTestsets()

# Map trainingsets and test sets to PyBrain
DS = ClassificationDataSet(trainingset['input_dimension'], trainingset['output_dimension'])
for i in range(0, len(trainingset['input_arrays'])):
	DS.appendLinked(trainingset['input_arrays'][i] , trainingset['output_arrays'][i])
DStest = ClassificationDataSet(trainingset['input_dimension'], trainingset['output_dimension'])
for i in range(0, len(testset['input_arrays'])):
	DStest.appendLinked(testset['input_arrays'][i], testset['output_arrays'][i])

# Create network
fnn = buildNetwork( DS.indim, 50, DS.outdim, outclass=SoftmaxLayer, fast=False )
#fnn = buildNetwork( DS.indim, 5, DS.outdim, outclass=SoftmaxLayer )
# Create trainer
trainer = BackpropTrainer(fnn, dataset=DS, momentum=0.01, verbose=True, weightdecay=0.0001)

#Train network and plot
if not args.l:
	errorlist = trainer.trainUntilConvergence(maxEpochs=args.e,validationProportion=0.2)
	if args.p:
		validation_plot(errorlist)
#Load network
if args.l:
        f = open (args.l, "r")
        fnn = pickle.load(f)

#Save network 
if args.s:
	with open("network.obj", "w") as f:
		pickle.dump(fnn, f)

out = fnn.activateOnDataset(DStest)

cat = collection.getUniqueCategories()
length = len(collection.getUniqueCategories())
cat = list(cat)
print "\n"+'Category'+"\t"+'Result'
for i in range(0, (len(collection.getUniqueCategories()) -1)):
#	print cat[i]+"\t"+str(out[0][i])
	print cat[i]
	for j in range(0,len(out)-1):
		print str(out[j][i])
#print cat
#pprint (out)
