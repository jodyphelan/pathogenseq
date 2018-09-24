#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
    bcf = ps.bcf(args.bcf_file)
    bcf.itol_from_bcf(args.mutation_file,args.amino_acid,args.no_ref,args.no_missing)

parser = argparse.ArgumentParser(description='bcf2matrix.py',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf_file', help='First read file')
parser.add_argument('mutation_file',default=None,type=str, help='First read file')
parser.add_argument('--amino_acid',action='store_true')
parser.add_argument('--no_ref',action='store_true')
parser.add_argument('--no_missing',action='store_true')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
