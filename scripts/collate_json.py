import sys
import json

sample_file = sys.argv[1]
extension = sys.argv[2]

samples = [l.rstrip() for l in open(sample_file).readlines()]
results = {}
for s in samples:
	results[s] = json.load(open("%s.%s" % (s,extension)))

keys = results[s].keys()

print "\t".join(keys)
for s in samples:
	print "\t".join([str(results[s][k]) for k in keys])
