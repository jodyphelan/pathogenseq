import pathogenseq as ps
import sys
import argparse

def main(args):
	fasta = ps.fasta(args.fasta)
	fasta.add_meta_data(args.data_file,args.outfile,args.delimiter)


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('fasta',help='bcf file')
parser.add_argument('data_file',help='reference file')
parser.add_argument('outfile',help='reference file')
parser.add_argument('--delimiter',default="_",type=str,help='reference file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
