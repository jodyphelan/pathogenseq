#! /usr/bin/env python
import pathogenseq as ps
import sys

if len(sys.argv)!=5:
	ps.log("sambamba_depth.py <ref> <bam> <prefix> <threads>",True)


ref_file = sys.argv[1]
bam_file = sys.argv[2]
prefix = sys.argv[3]
threads = sys.argv[4]


bam = ps.bam(bam_file,prefix,ref_file,threads=threads)
out_file = "%s.sambamba.depth" % prefix
bam.sambamba_depth(out_file)
