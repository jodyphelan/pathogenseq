#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
	bcf_file = args.bcf
	ref_file = args.ref
	bcf = ps.bcf(bcf_file)
	bcf.generate_consensus(ref_file)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf',help='bcf file')
parser.add_argument('ref',help='reference file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
