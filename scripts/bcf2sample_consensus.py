import pathogenseq as ps
import sys



import argparse

def main(args):
	bcf = ps.bcf(args.vcf)
	bcf.per_sample_bcf2fa(args.sample,args.ref,args.no_chrom)


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('vcf',help='bcf file')
parser.add_argument('sample',help='bcf file')
parser.add_argument('ref',help='reference file')
parser.add_argument('--no-chrom',action="store_true",help='reference file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
