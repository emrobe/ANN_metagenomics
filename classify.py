from pybrain.datasets            import ClassificationDataSet
from pybrain.utilities           import percentError
from pybrain.tools.shortcuts     import buildNetwork
from pybrain.tools.validation    import CrossValidator
from pybrain.supervised.trainers import BackpropTrainer
from pybrain.structure.modules   import SoftmaxLayer
from pybrain.tools.validation    import ModuleValidator

from pylab import ion, ioff, figure, draw, contourf, clf, show, hold, plot
from scipy import diag, arange, meshgrid, where
from numpy.random import multivariate_normal

from biom_parser import biom_parser as bp
from collection import collection as col
from plot_validation import validation_plot
from plot_testsets import testset_plot 
import os, sys
from pprint import pprint

#import arac
import argparse
import pickle 

parser = argparse.ArgumentParser(description=\
  'This program classifies taxonomic classifications with a PyBrian neural network')
parser.add_argument('-e',type=int, required=False,\
  help="Number of epochs to train")
#parser.add_argument('-pv',action="store_true", default=False, required=False,\
#  help="Plot training and validation error")
parser.add_argument('-pt',action="store_true", default=False, required=False,\
  help="Plot heatmap of testsets")
parser.add_argument('-s',action="store_true",required=False,\
  help="Use pickle to save network to file for future use")
parser.add_argument('-l',type=str,required=False,\
  help="Load a trained network saved with pickle")
parser.add_argument('--train',type=str,required=False,\
  help="Input firectory for trainingsets")
parser.add_argument('--test',type=str,required=False,\
  help="Input directory for test sets")
args = parser.parse_args()


#traindir = "train"
#testdir = "test"
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
for file in os.listdir(args.test):
    if file.endswith('fasta'):
		a = bp(os.path.join(args.test,file), os.path.join(args.test, 'metadata.csv'))
		a.createMgrast()
		a.parseMgrastCategories('3', '1')
        #Level 1 = Root, Level 2 = Kingdom, Level 3 = Phylum...etc
		#a.set_lowest_taxa_level(lower_taxa_level)
	        #a.delete_upper_taxa_level(upper_taxa_level)
	        #a.createDataDict()
		#a.convertToRelative()
		#a.convertToShares()
		a.convertToPresence()
		for otu in a.getOtuList():
			collection.addTaxa(otu)
		collection.addFile(a, 'test')

# Calculate and print number of total inputnodes (unique taxa) and total output nodes (uniqe categories to classify)
collection.setUniqueTaxa()
collection.setUniqueCategories()
print 'Unique Taxa (#input nodes): '+str(len(collection.getUniqueTaxa()))
print 'Unique Categories (#output nodes): '+str(len(collection.getUniqueCategories()))

# Create trainingsets and test sets
trainingset = collection.createAnnTrainingsets()
	

# Map trainingsets and test sets to PyBrain
DS = ClassificationDataSet(trainingset['input_dimension'], trainingset['output_dimension'])
for i in range(0, len(trainingset['input_arrays'])):
	DS.appendLinked(trainingset['input_arrays'][i] , trainingset['output_arrays'][i])


fnn = buildNetwork( DS.indim, 100, DS.outdim, outclass=SoftmaxLayer, fast=False )

trainer = BackpropTrainer(fnn, dataset=DS, momentum=0.01, verbose=True, weightdecay=0.0001)

#Train network and plot
if not args.l:
	errorlist = trainer.trainEpochs(args.e)

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

#Plot heatmap
if args.pt:
        testset_plot(xcategories, ytestsets, zdata)

