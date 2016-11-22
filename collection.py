class collection():

	def __init__(self):
		self.test = {}
		self.testAnnArray = {'in':[], 'out':[]}
		self.train = {}
		self.trainAnnArray = {'in':[], 'out':[]}
		self.all_taxa = []
		self.unique_taxa = []
		self.all_categories = []
		self.unique_categories = []

	def addFile(self, a, filetype):
		if filetype == 'train':
			self.train[a.getName()] = a
		elif filetype == 'test':
			self.test[a.getName()] = a
		else:
			exit("Wrong filetype "+filetype)

	def addTaxa(self, string):
		self.all_taxa.append(string)

	def addCategory(self, category):
		self.all_categories.append(category)

	def setUniqueTaxa(self):
		self.unique_taxa = set(self.all_taxa)

	def getUniqueTaxa(self):
		return self.unique_taxa

	def setUniqueCategories(self):
		self.unique_categories = set(self.all_categories)

	def getUniqueCategories(self):
		return self.unique_categories

	def printtmp(self):
		print self.unique_categories

	def createAnnTrainingsets(self):
		ds = {}
		ds['input_dimension'] = len(self.unique_taxa)
		ds['output_dimension'] = len(self.unique_categories)
		ds['input_arrays'] = []
		ds['output_arrays'] = []
		for mastertaxa in self.unique_taxa:
			tmparray = []
			for biom_object in self.train:
				if mastertaxa in self.train[biom_object].taxa_counts:
					tmparray.append(self.train[biom_object].taxa_counts[mastertaxa])
				else:
					tmparray.append(0)
			ds['input_arrays'].append(tmparray)
		for mastercategory in self.unique_categories:
			tmparray = []
			for biom_object in self.train:
				if mastercategory in self.train[biom_object].categories[0]:
					tmparray.append(1)
				else:
					tmparray.append(0)
			ds['output_arrays'].append(tmparray)
		transposedinput = map(list, zip(*ds['input_arrays']))
		transposedoutput = map(list, zip(*ds['output_arrays']))
		ds['input_arrays'] = transposedinput
		ds['output_arrays'] = transposedoutput
		#print ds['input_dimension']
		#print ds['output_dimension']
		#print ds['input_arrays']
		#print ds['output_arrays']
		return ds

	def createAnnTestsets(self):
                ds = {}
                ds['input_dimension'] = len(self.unique_taxa)
                ds['output_dimension'] = len(self.unique_categories)
                ds['input_arrays'] = []
                ds['output_arrays'] = []
                for mastertaxa in self.unique_taxa:
                        tmparray = []
                        for biom_object in self.test:
                                if mastertaxa in self.test[biom_object].taxa_counts:
                                        tmparray.append(self.test[biom_object].taxa_counts[mastertaxa])
                                else:
                                        tmparray.append(0)
                        ds['input_arrays'].append(tmparray)
                for mastercategory in self.unique_categories:
                        tmparray = []
                        for biom_object in self.test:
                                if mastercategory in self.test[biom_object].categories[0]:
                                        tmparray.append(1)
                                else:
                                        tmparray.append(0)
                        ds['output_arrays'].append(tmparray)
                transposedinput = map(list, zip(*ds['input_arrays']))
                transposedoutput = map(list, zip(*ds['output_arrays']))
                ds['input_arrays'] = transposedinput
                ds['output_arrays'] = transposedoutput
                #print ds['input_dimension']
                #print ds['output_dimension']
                #print ds['input_arrays']
                #print ds['output_arrays']
                return ds

	def createAnnTestset(self,testset):
                ds = {}
                ds['input_dimension'] = len(self.unique_taxa)
                ds['output_dimension'] = len(self.unique_categories)
                ds['input_arrays'] = []
                ds['output_arrays'] = []
                for mastertaxa in self.unique_taxa:
                        tmparray = []
			if mastertaxa in self.test[testset].taxa_counts:
				tmparray.append(self.test[testset].taxa_counts[mastertaxa])
			else:
				tmparray.append(0)
                        ds['input_arrays'].append(tmparray)
                for mastercategory in self.unique_categories:
                        tmparray = []
			if mastercategory in self.test[testset].categories[0]:
				tmparray.append(1)
			else:
				tmparray.append(0)
                        ds['output_arrays'].append(tmparray)
		#print ds['input_arrays']
                transposedinput = map(list, zip(*ds['input_arrays']))
                transposedoutput = map(list, zip(*ds['output_arrays']))
                ds['input_arrays'] = transposedinput
                ds['output_arrays'] = transposedoutput
                #print ds['input_dimension']
                #print ds['output_dimension']
                #print ds['input_arrays']
                print ds['output_arrays']
                return ds
