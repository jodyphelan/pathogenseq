#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse
import json

def main(args):
	bam = ps.bam(args.bam,args.prefix,args.ref,threads=args.threads)
	bcf = bam.call_variants(call_method=args.call_method,gff_file=args.gff,bed_file=args.bed,mixed_as_missing=False,threads=args.threads,min_dp=args.min_depth)
	csq = bcf.load_csq_alt(ann_file=args.ann,changes=True)
	tmp_bcf = "%s.missing.bcf" % args.prefix
	missing_pos = ps.get_missing_positions(tmp_bcf)
	outfile = "%s.results.json" % args.prefix
	results = {"variants":[],"missing":missing_pos}
	for sample in csq:
		results["variants"]  = csq[sample]
	if args.barcode:
		mutations = bam.get_bed_gt(args.barcode)
		barcode = ps.barcode(mutations,"lineages.bed")
		results["barcode"] = barcode


	json.dump(results,open(outfile,"w"))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bam', help='BAM file')
parser.add_argument('ref', help='Reference Sequence')
parser.add_argument('gff', help='GFF file')
parser.add_argument('bed', help='BED file')
parser.add_argument('prefix', help='Prefix for files')
parser.add_argument('--ann','-a',type=str,default=None, help='Annotation file')
parser.add_argument('--barcode','-b',type=str,default=None, help='Barcode bed file')
parser.add_argument('--threads','-t', type=int,default=1,help='Number of threads')
parser.add_argument('--call_method','-c', type=str,default="low",choices=["low","high","optimise"],help='Number of threads')
parser.add_argument('--min_depth','-d', type=int,default=10,help='Number of threads')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
