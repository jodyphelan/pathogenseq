#! /usr/bin/env python
import pathogenseq as ps
import argparse
import json

def main(args):
    if not args.r1:
        ps.log("Please provide at least one fastq file with -1...Exiting")
        quit()
    else:
        ps.filecheck(args.r1)
    if args.r2: ps.filecheck(args.r2)
    if not args.prefix:
        ps.log("Please provide a file output prefix...Exiting")
        quit()

    stats = {}
    fastqqc = ps.qc_fastq(args.prefix,args.r1,args.r2) if args.r2 else ps.qc_fastq(args.prefix,args.r1)
    stats["mean_read_len"] = fastqqc.mean_read_len
    stats["read_num"] = fastqqc.read_num

    stats_json = "%s.fastq_stats.json" % args.prefix
    json.dump(stats,open(stats_json,"w"))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--r1','-1', help='First read file')
parser.add_argument('--r2','-2', help='Second read file')
parser.add_argument('--prefix','-p', help='Prefix for files')
parser.add_argument('--kraken', help='Use kraken files')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
