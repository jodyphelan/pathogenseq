import pathogenseq as ps
import sys
import argparse

def main(args):
	tree = ps.tree(args.tree)
	tree.rename_nodes(args.index_file,args.outfile,args.strict,args.append)


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('tree',help='bcf file')
parser.add_argument('index_file',help='reference file')
parser.add_argument('outfile',help='reference file')
parser.add_argument('--strict',action="store_true",help='reference file')
parser.add_argument('--append',default=None,type=str,help='reference file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
