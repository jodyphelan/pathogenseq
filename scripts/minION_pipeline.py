#! /usr/bin/env python
import sys
import pathogenseq as ps
import json
from collections import OrderedDict

ref = sys.argv[1]
r1 = sys.argv[2]
prefix = sys.argv[3]
threads = sys.argv[4]

cov_png = "%s.cov.png" % prefix
stats_file = "%s.stats.json" % prefix
gc_file = "%s.gc_skew.txt" % prefix
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
bam_qc.plot_cov("Chromosome",cov_png)
bam_qc.extract_gc_skew(gc_file)
variants = bam.pileup2vcf()
stats["hom_variants"] = len([x for x in variants if x[5]=="1/1"])
stats["het_variants"] = len([x for x in variants if x[5]=="0/1"])
json.dump(stats,open(stats_file,"w"))
