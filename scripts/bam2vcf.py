#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
	bam = ps.bam(bam_file=args.bam,ref_file=args.ref,prefix=args.prefix,threads=args.threads)
	bam.gbcf(threads=args.threads,platform=args.platform)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bam', help='First read file')
parser.add_argument('ref', help='Second read file')
parser.add_argument('prefix', help='Reference Sequence')
parser.add_argument('--threads','-t', type=int, default=1, help='Number of threads')
parser.add_argument('--platform','-p', type=str,default="illumina",choices=["illumina","minION"],help='Mapping tool to use')
parser.add_argument('--primers', type=str,default=None,help='Mapping tool to use')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
