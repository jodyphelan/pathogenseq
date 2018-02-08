import sys
import json

j = json.load(open(sys.argv[1]))
columns = sorted(j[j.keys()[0]].keys())
print "sample\t%s" % "\t".join(columns)
for s in j:
	print "%s\t%s" % (s,"\t".join([j[s][c] for c in columns]))
