#! /usr/bin/env python
from __future__ import division
import sys
import pathogenseq as ps
import argparse
import json

def main(args):
	ref = args.ref
	r1 = args.r1
	r2 = args.r2
	prefix = args.prefix
	threads = args.threads

	stats = {}
	fastqqc = ps.qc_fastq(prefix,r1,r2)
	stats["fastq_mean_read_len"] = fastqqc.mean_read_len
	stats["fastq_read_num"] = fastqqc.read_num


	fastq = ps.fastq(prefix,ref,r1,r2,threads=threads)

	fastq.illumina(mapper="minimap2")

	bam_file = "%s.bam" % prefix

	bam = ps.bam(bam_file,prefix,ref)
	bam.gbcf(vtype="snps",threads=threads,call_method=args.call_method)

	bamqc = bam.get_bam_qc()
	cov_plot = "%s.cov.png" % (prefix)
	bamqc.plot_cov("Chromosome",cov_plot)
	stats["bam_pct_reads_mapped"] = bamqc.pct_reads_mapped
	stats["bam_med_dp"] = bamqc.med_dp
	stats["bam_depth_10"] = bamqc.genome_cov[10]
	stats["bam_depth_5"] = bamqc.genome_cov[5]

	json.dump(stats,open("%s.stats.json" % prefix,"w"))
	O = open("%s.log"%prefix,"w")
	for x in ["fastq_mean_read_len","fastq_read_num","bam_pct_reads_mapped","bam_med_dp","bam_depth_5","bam_depth_10"]:
		O.write("%s\t%s\n" % (x,stats[x]))
	O.close()


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('ref', help='First read file')
parser.add_argument('r1', help='First read file')
parser.add_argument('r2', help='First read file')
parser.add_argument('prefix', help='First read file')
parser.add_argument('--threads',"-t",type=int,default=1, help='First read file')
parser.add_argument('--call_method',"-b",type=str,default="optimise",choices=["optimise","low","high"],help='First read file')
parser.add_argument('--centrifuge',"-c",default=None, help='First read file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
