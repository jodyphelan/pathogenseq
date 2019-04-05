#! /usr/bin/env python
import pathogenseq as ps
import sys
import argparse

def main(args):
	ps.mauve_call_variants(args.ref,args.query,args.prefix)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('ref',help='NGS Platform')
parser.add_argument('query',help='NGS Platform')
parser.add_argument('prefix',help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
