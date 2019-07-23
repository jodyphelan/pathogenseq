#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
	bcf = ps.bcf(args.vcf)
	bcf.vcf_to_fasta(outfile=args.outfile,ref_file=args.ref,threads=args.threads)


parser = argparse.ArgumentParser(description='bcf2fasta.py',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('vcf', help='VCF file')
parser.add_argument('ref', type=str, help='Reference')
parser.add_argument('outfile', type=str, help='Output name')
parser.add_argument('--threads', type=int, default=4, help='Number of threads')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
