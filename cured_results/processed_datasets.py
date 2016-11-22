import sys
table = open('../table.tsv', 'r')
for line in table:
	id = line.strip().split()
	if id[2] in sys.argv:
		print line.strip()
