#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
	bcf = ps.bcf(args.vcf)
	bcf.distance(args.outfile)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf','-v', help='VCF file',required=True)
parser.add_argument('--outfile','-o', help='OUput file',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
