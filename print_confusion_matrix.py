import pickle
import pandas as pd

f = open ('confusion_matrix_e2000n30.obj', "r")
conf = pickle.load(f)
print dir(conf)
for i in conf.iteritems():
	print i
with open ('confusion_matrix_e2000n30.csv', "w") as csv:
	print conf.to_csv(csv)
