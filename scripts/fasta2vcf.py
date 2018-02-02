#! /usr/bin/env python
import sys
import pathogenseq as ps

infile = sys.argv[1]
ref_file = sys.argv[2]
prefix = sys.argv[3]
fasta = ps.fasta(infile)
fasta.get_VCF(ref_file,prefix)
