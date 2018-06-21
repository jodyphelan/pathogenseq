#! /usr/bin/env python
import sys
import json

sample_file = sys.argv[1]
extension = sys.argv[2]

samples = [l.rstrip() for l in open(sample_file).readlines()]
results = {}
keys = set()
for s in samples:
	results[s] = json.load(open("%s%s" % (s,extension)))
	for x in results[s]:
		keys.add(x)

print("sample\t%s"%"\t".join(keys))
for s in samples:
	print("%s\t%s" % (s,"\t".join([str(results[s][k]) if k in results[s] else "-" for k in keys])))
