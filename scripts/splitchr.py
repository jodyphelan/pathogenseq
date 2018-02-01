#! /usr/bin/env python
import sys
import pathogenseq as ps

fasta = ps.fasta(sys.argv[1])
fasta.splitchr(int(sys.argv[2]))
