#! /usr/bin/env python
import pathogenseq
import argparse

def main(args):
	if not args.prefix:
		pathogenseq.log("Please specify prefix with -p")
		quit(1)
	if not args.ref:
		log("Please use --ref to provide a reference... Exiting",ext=T)
	x= pathogenseq.fastq(args.prefix,args.ref,args.r1,args.r2,threads=args.threads)
	x.illumina(mapper=args.mapper)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--r1','-1', help='First read file')
parser.add_argument('--r2','-2', help='Second read file')
parser.add_argument('--ref','-r', help='Reference Sequence')
parser.add_argument('--threads','-t', type=int, default=1, help='Number of threads')
parser.add_argument('--prefix','-p', help='Prefix for files')
parser.add_argument('--mapper','-m', type=str,default="bwa",choices=["bwa","minimap2","bowtie2"],help='Mapping tool to use')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
