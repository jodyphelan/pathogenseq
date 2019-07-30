#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
	bcf_file = args.bcf
	ref_file = args.ref
	bcf = ps.bcf(bcf_file)
	bcf.generate_consensus(ref_file,threads=args.threads,no_chrom=args.combine)
	if args.combine:
		files = " ".join(["%s.%s.fasta" % (bcf.prefix,s) for s in bcf.samples])
		ps.run_cmd("cat %s > %s.genome.fasta" % (files,bcf.prefix))
		ps.run_cmd("rm %s" % (files))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf',help='bcf file')
parser.add_argument('ref',help='reference file')
parser.add_argument('--threads','-t',default=4,type=int,help='reference file')
parser.add_argument('--combine',action="store_true",help='reference file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
