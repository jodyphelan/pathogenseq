#! /usr/bin/env python
import sys
import pathogenseq as ps

if len(sys.argv)!=4:
	print("examl.py <snps.fasta> <prefix> <threads>")
	quit()

snps_file = sys.argv[1]
prefix = sys.argv[2]
threads = sys.argv[3]
phylo = ps.phylo(snps_file,prefix,threads)
phylo.examl()
