#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse
import json

def main(args):
	bam = ps.bam(args.bam,args.prefix,args.ref)
	bcf = bam.call_variants(call_method=args.depth,gff_file=args.gff,bed_file=args.bed,mixed_as_missing=True,threads=args.threads)
	csq = bcf.load_csq(ann_file=args.ann,changes=True)
	tmp_bcf = "%s.missing.bcf" % args.prefix
	missing_pos = ps.get_missing_positions(tmp_bcf)
	outfile = "%s.results.json" % args.prefix
	results = {"variants":{},"missing":missing_pos}
	for gene in csq:
		results["variants"][gene] = []
		for var in csq[gene]:
			results["variants"][gene].append(var.values()[0])
	json.dump(results,open(outfile,"w"))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bam', help='BAM file')
parser.add_argument('ref', help='Reference Sequence')
parser.add_argument('gff', help='GFF file')
parser.add_argument('bed', help='BED file')
parser.add_argument('ann', help='ANN file')
parser.add_argument('prefix', help='Prefix for files')
parser.add_argument('--threads','-t', type=int,default=1,help='Number of threads')
parser.add_argument('--depth','-d', type=str,default="optimise",choices=["low","high","optimise"],help='Number of threads')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
