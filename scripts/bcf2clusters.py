#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
    bcf = ps.bcf(args.bcf_file)
    clusters = bcf.get_clusters(args.dist,args.meta,args.cols,args.remove_singletons)


parser = argparse.ArgumentParser(description='bcf2matrix.py',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf_file', help='First read file')
parser.add_argument('--dist','-d', default=10,type=int)
parser.add_argument('--meta','-m', default=None,type=str)
parser.add_argument('--cols','-c', default=None,type=str)
parser.add_argument('--remove_singletons','-r', action="store_true")
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
