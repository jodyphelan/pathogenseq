#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):

	x = ps.qc_bam(args.bam,args.ref)
	x.plot_cov(args.chrom,args.outfile,start=args.start,end=args.end,window=args.window,step=args.step,optimise=args.optimise)



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bam', help='BAM file')
parser.add_argument('ref', help='Reference file')
parser.add_argument('chrom', help='Chromosome')
parser.add_argument('outfile', help='Chromosome')
parser.add_argument('--start','-s', type=int, default=None, help='Number of threads')
parser.add_argument('--end','-e', type=int, default=None, help='Prefix for files')
parser.add_argument('--window','-w', type=int, default=10000, help='Prefix for files')
parser.add_argument('--step','-x', type=int, default=5000, help='Prefix for files')
parser.add_argument('--optimise', action="store_true", help='Prefix for files')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
