#! /usr/bin/env python
import sys
import pathogenseq as ps
import json
from collections import OrderedDict
import argparse

def main(args):
	ref = args.ref
	r1 = args.reads
	prefix = args.prefix
	threads = args.threads

	cov_png = "%s.cov.png" % prefix
	stats_file = "%s.stats.json" % prefix
	gc_file = "%s.gc_skew.json" % prefix
	cov_file = "%s.regions.cov.json" % prefix
	stats = OrderedDict()
	fq = ps.fastq(prefix,ref,r1,threads=threads)
	fq_qc = fq.get_fastq_qc()
	stats["mean_read_len"] = fq_qc.mean_read_len
	stats["median_read_len"] = fq_qc.median_read_len
	stats["read_num"] = fq_qc.read_num
	bam = fq.minION()
	bam_qc = bam.get_bam_qc()
	stats["med_dp"] = bam_qc.med_dp
	stats["pct_reads_mapped"] = bam_qc.pct_reads_mapped
	stats["genome_cov_1"] = bam_qc.genome_cov[1]
	stats["genome_cov_10"] = bam_qc.genome_cov[10]
	fasta = ps.fasta(ref).fa_dict
	for seq in fasta:
		bam_qc.plot_cov(seq,cov_png,primers=args.primers)
	bam_qc.extract_gc_skew(gc_file)
	if args.bed_cov: bam_qc.save_cov(cov_file,args.bed_cov)
	variants = bam.pileup2vcf(indels=False)
	stats["hom_variants"] = len([x for x in variants if x[5]=="1/1"])
	stats["het_variants"] = len([x for x in variants if x[5]=="0/1"])
	json.dump(stats,open(stats_file,"w"))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('ref', help='First read file')
parser.add_argument('reads', help='First read file')
parser.add_argument('prefix', help='First read file')
parser.add_argument('--threads',"-t", help='First read file')
parser.add_argument('--bed_cov',"-b",default=None, help='First read file')
parser.add_argument('--primers',"-p",default=None, help='First read file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
