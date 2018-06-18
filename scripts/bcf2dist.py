#! /usr/bin/env python
import sys
import pathogenseq as ps

infile,outfile = sys.argv[1:]

bcf = ps.bcf(infile)
bcf.distance(outfile)
