#import biom
import collections
import os

class biom_parser():

	def __init__(self, inputfile, categories_file):
		self.inputfile = inputfile
		self.categories_file = categories_file
		self.tax_level_list = []
		self.taxa_counts = collections.defaultdict(int)
		self.taxa_share_level = collections.defaultdict(float)
		self.categories = []

	def parseCategories(self):
		categories_file = open(self.categories_file)
		filename = os.path.split(self.inputfile)
		for line in categories_file:
			tmp = line.strip().split("\t")
			if tmp[0] == filename[-1]:
				tmp.pop(0)
				self.categories.append(tmp)
		if not self.categories:
			exit(filename[-1]+" did not parse any categories!")

	def parseMgrastCategories(self, columns, id_column):
		categories_file = open(self.categories_file)
		filename = os.path.split(self.inputfile)

		for line in categories_file:
			tmp = line.strip().split("\t")
			if tmp[int(id_column)] == filename[-1]:
				categories = []
				for column in columns:
					categories.append(tmp[int(column)])
				self.categories.append(categories)
		if not self.categories:
			exit(filename[-1]+" did not parse any categories!")

	def convertToPresence(self):
		for taxa in self.taxa_counts:
			self.taxa_counts[taxa] = 1

	def convertToShares(self):
		for taxa in self.taxa_counts:
			self.taxa_counts[taxa] = self.taxa_share_level[taxa]

	def convertToRelative(self):
		total = 0
		for taxa in self.taxa_counts:
			total += int(self.taxa_counts[taxa])
		for taxa in self.taxa_counts:
			self.taxa_counts[taxa] = float(self.taxa_counts[taxa]) / float(total)

	def createDataDict(self):
		if not self.tax_level_list:
			exit("No taxa in list! createEbiBiomList() first.")
		for otu in self.tax_level_list:
			for taxa in otu:
				self.taxa_counts[taxa] += 1

	def getName(self):
		return self.inputfile

	def createEbiBiomList(self):
		tax_level_list = []
		table = biom.load_table(self.inputfile)
		length = len(table.ids(axis='observation'))
		#print table.metadata('0', axis='observation').taxonomy()
		for i in range (0, length):
			id = str(i)
			a = table.metadata(id, axis='observation')
			self.tax_level_list.append(a['taxonomy'])

	def createMgrast(self):
		#This method is mutually exclusive to createDataDict() and createEbiBiomList()!
		tax_level_list = []
		data = open(self.inputfile)
		for line in data:
		    line = line.split('\t')
		    if (len(line) == 6):
		        if line[1] == 'Taxa':
		            continue
		        self.tax_level_list.append(line[1])
		        self.taxa_counts[line[1]] = line[2]
			self.taxa_share_level[line[1]] = line[3]

	def delete_upper_taxa_level(self, level):
		tmp = self.tax_level_list
		for otu in tmp:
			for i in range(1,level):
				if len(otu) != 0:
					otu.pop(0)
		tmp2 = [x for x in tmp if x != []]
		self.tax_level_list = tmp2

        def set_lowest_taxa_level(self, level):
                tmp = self.tax_level_list
                for otu in tmp:
			count = (len(otu)-level)
                        if count>=1:
				for i in range(0,count):
					otu.pop()
                self.tax_level_list = tmp

	def printList(self):
		print self.tax_level_list

	def printAll(self):
		print self.taxa_counts
		print self.tax_level_list
		print self.taxa_share_level

	def getOtuList(self):
		return self.tax_level_list

	def getCategoryList(self):
		return self.categories

	def selectTaxaRange(self, upper, lower):
                _set_lowest_taxa_level(self, lower)
                _delete_upper_taxa_level(self, upper)
