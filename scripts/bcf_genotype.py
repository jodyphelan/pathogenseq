#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
	bcf = ps.bcf(args.bcf)
	print(bcf.get_bed_gt(args.bed,args.ref))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf',help='bcf file')
parser.add_argument('bed',help='bcf file')
parser.add_argument('ref',help='reference file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
