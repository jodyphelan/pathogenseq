#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse

def main(args):
    bcf = ps.bcf(args.bcf_file)
    bcf.get_plink_dist()

parser = argparse.ArgumentParser(description='bcf2matrix.py',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf_file', help='First read file')
parser.add_argument('--no_missing',action='store_true')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
