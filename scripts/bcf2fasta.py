#! /usr/bin/env python
import sys
import pathogenseq as ps
infile = sys.argv[1]
reffile = sys.argv[2]
outfile = sys.argv[3]
threads = sys.argv[4]

bcf = ps.bcf(infile)
bcf.vcf_to_fasta_alt(outfile=outfile,ref_file=reffile,threads=threads)
