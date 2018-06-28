#! /usr/bin/env python
import pathogenseq as ps
import argparse
import json

def main(args):
    if not args.bam:
        ps.log("Please provide bam file")
        quit()
    else:
         ps.filecheck(args.bam)
    if not args.ref:
        ps.log("Please provide reference")
        quit()
    else:
        ps.filecheck(args.ref)
        fasta = ps.fasta(args.ref)
    if not args.prefix:
        ps.log("Please provide prefix")
        quit()
    if args.gff and args.bed:
        ps.log("Please provide either a GFF file or BED file but not both...Exiting!")
        quit()
    if args.gff:
        if not args.gffkey:
            ps.log("Please provide the key to look for in the GFF file...Exiting!")
            quit()
        else:
            ps.filecheck(args.gff)
    if args.bed: ps.filecheck(args.bed)

    cov_json = "%s.cov.json" % args.prefix
    stats_json = "%s.bam_stats.json" % args.prefix
    region_json = "%s.regions.cov.json" % args.prefix

    stats = {}

    bamqc = ps.qc_bam(args.bam,args.ref)
    for s in fasta.fa_dict:
        cov_plot = "%s.%s.cov.png" % (args.prefix,s)
        bamqc.plot_cov(s,cov_plot)
    bamqc.save_cov(cov_json)
    stats["pct_reads_mapped"] = bamqc.pct_reads_mapped
    stats["med_dp"] = bamqc.med_dp
    region_cov = {}
    if args.gff:
        region_cov = bamqc.gff_cov(args.gff,args.gffkey)
    elif args.bed:
        region_cov = bamqc.bed_cov(args.bed)
    json.dump(stats,open(stats_json,"w"))
    json.dump(region_cov,open(region_json,"w"))





parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam','-a', help='BAM file')
parser.add_argument('--ref','-r', help='Reference Sequence')
parser.add_argument('--prefix','-p', help='Prefix for files')
parser.add_argument('--gff','-g', help='GFF file')
parser.add_argument('--gffkey','-k', help='GFF key')
parser.add_argument('--bed','-b', help='BED file')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
