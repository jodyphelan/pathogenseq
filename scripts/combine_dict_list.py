#! /usr/bin/env python
import sys
import json

prefix_file = sys.argv[1]
extension = sys.argv[2]
prefixes = [l.rstrip() for l in open(prefix_file).readlines()]

for pre in prefixes:
	f = "%s%s" % (pre,extension)
	j = json.load(open(f))
	for d in j:
		for x in j[d]:
			print("%s\t%s\t%s" % (pre,d,x))
