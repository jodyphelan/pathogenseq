#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
    bcf = ps.bcf(args.bcf_file)
    bcf.extract_matrix(args.out,args.format,args.annotation)

parser = argparse.ArgumentParser(description='bcf2matrix.py',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf_file', help='First read file')
parser.add_argument('--out','-o',default=None,type=str, help='First read file')
parser.add_argument('--format','-f',default="old",choices=["old","new"], type=str, help='First read file')
parser.add_argument('--annotation', action="store_true", help='First read file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
