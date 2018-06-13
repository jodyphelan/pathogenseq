#! /usr/bin/env python
import sys
import pathogenseq as ps
import json
from collections import OrderedDict
import argparse

def main(args):
	ref = args.ref
	r1 = args.r1
	r2 = args.r2
	prefix = args.prefix
	threads = args.threads

	cov_png = "%s.cov.png" % prefix
	stats_file = "%s.stats.json" % prefix
	gc_file = "%s.gc_skew.json" % prefix
	cov_file = "%s.regions.cov.json" % prefix
	stats = OrderedDict()
	fq = ps.fastq(prefix,ref,r1,r2,threads=threads)
	fq_qc = fq.get_fastq_qc()
	if args.centrifuge:
		t1,t2 = fq_qc.run_centrifuge(args.centrifuge,False,threads)
		stats["centrifuge_top_hit"] = t1
		stats["centrifuge_top_hit_num_reads"] = t2
	stats["mean_read_len"] = fq_qc.mean_read_len
	stats["median_read_len"] = fq_qc.median_read_len
	stats["read_num"] = fq_qc.read_num
	bam = fq.illumina(mapper=args.mapper)
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
	variants = bam.gbcf(primers=args.primers,chunk_size=args.window)
	bcfstats = variants.load_stats()
	stats["hom_variants"] = bcfstats["PSC"][prefix]["nNonRefHom"]
	stats["het_variants"] = bcfstats["PSC"][prefix]["nHets"]
	stats["hom_ref"] = bcfstats["PSC"][prefix]["nRefHom"]
	json.dump(stats,open(stats_file,"w"))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('ref', help='First read file')
parser.add_argument('r1', help='First read file')
parser.add_argument('r2', help='First read file')
parser.add_argument('prefix', help='First read file')
parser.add_argument('--threads',"-t",type=int,default=1, help='First read file')
parser.add_argument('--mapper',"-m",type=str,choices=["bwa","minimap2","bowtie2"],default="bwa", help='First read file')
parser.add_argument('--bed_cov',"-b",default=None, help='First read file')
parser.add_argument('--primers',"-p",default=None, help='First read file')
parser.add_argument('--centrifuge',"-c",default=None, help='First read file')
parser.add_argument('--window',default=50000,type=int, help='First read file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
