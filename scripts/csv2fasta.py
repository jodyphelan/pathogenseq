#! /usr/bin/env python
import csv
import sys
infile = sys.argv[1]
outfile = sys.argv[2]
OUT = open(outfile,"w")
for row in csv.DictReader(open(infile)):
	OUT.write(">%(ID)s\n%(Sequence)s\n" % row)
OUT.close()
