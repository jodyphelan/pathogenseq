import sys
import csv
from collections import defaultdict

meta = defaultdict(dict)
samples = []
for row in csv.DictReader(open(sys.argv[1])):
	print row
	meta_columns = set(row.keys())-set(["sample"])
	for x in meta_columns:
		meta[row["sample"]][x] = row[x]
	samples.append(row["sample"])


for c in meta_columns:
	data = [meta[s][c] for s in samples]
	binary =  True if set(data) == set(["0","1",""]) or set(data)== set(["0","1"]) else False
	if binary:
		colour_dict = {"1":"#000000","0":"#ffffff","":"#848484"}
		for s in samples:
			print "%s\t%s" % (s,colour_dict[meta[s][c]])
