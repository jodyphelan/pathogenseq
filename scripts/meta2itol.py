import sys
import csv
from collections import defaultdict
from colour import Color

meta = defaultdict(dict)
samples = []
for row in csv.DictReader(open(sys.argv[1])):
	meta_columns = set(row.keys())-set(["sample"])
	for x in meta_columns:
		meta[row["sample"]][x] = row[x] if row[x] != "" else "N/A"
	samples.append(row["sample"])


for c in meta_columns:
	O = open("%s.meta.itol.txt" % c,"w")
	header = """DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Lineage
COLOR	#ff0000

LEGEND_TITLE	%s
""" % c
	O.write(header)


	raw_data = [meta[s][c] for s in samples]
	data = set(raw_data)-set(["N/A"])
	binary =  True if  data == set(["0","1"]) else False
	if binary:
		colour_dict = {"1":"#000000","0":"#ffffff","N/A":"#848484"}
	else:
		colours = list(Color("red").range_to(Color("blue"),len(data)))
		colour_dict = {d:colours[list(data).index(d)].get_hex() for d in data}
		colour_dict["N/A"]= "#848484"
	O.write("LEGEND_SHAPES\t%s\n" % "\t".join(["1" for _ in sorted(colour_dict)]))
	O.write("LEGEND_COLORS\t%s\n" % "\t".join([x for x in sorted(colour_dict)]))
	O.write("LEGEND_LABELS\t%s\n" % "\t".join([colour_dict[x] for x in sorted(colour_dict)]))
	O.write("\nDATA\n")
	for s in samples:
		O.write("%s\t%s\n" % (s,colour_dict[meta[s][c]]))
	O.close()
