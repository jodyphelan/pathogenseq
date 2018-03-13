#! /usr/bin/env python
from __future__ import division
import sys
import pathogenseq as ps
import argparse
import json

ref = sys.argv[1]
r1 = sys.argv[2]
r2 = sys.argv[3]
prefix = sys.argv[4]
threads = sys.argv[5]

stats = {}
fastqqc = ps.qc_fastq(prefix,r1,r2)
stats["fastq_mean_read_len"] = fastqqc.mean_read_len
stats["fastq_read_num"] = fastqqc.read_num

fastqqc.run_kraken("/opt/storage2/ernest/kraken/kraken/standard_db","77643,1773,78331,33894,1765")

fr1 = "%s_1.kraken_filt.fastq.gz" %prefix
fr2 = "%s_2.kraken_filt.fastq.gz" %prefix

newfastqqc = ps.qc_fastq(prefix,fr1,fr2)
stats["kraken_pct_filt"] = newfastqqc.read_num/fastqqc.read_num*100

fastq = ps.fastq(prefix,ref,fr1,fr2,threads=threads)

fastq.illumina(mapper="bowtie2")

bam_file = "%s.bam" % prefix

bam = ps.bam(bam_file,prefix,ref)
bam.gbcf(vtype="snps",threads=threads)

bamqc = bam.get_bam_qc(ref)
cov_plot = "%s.cov.png" % (prefix)
bamqc.plot_cov("Chromosome",cov_plot)
stats["bam_pct_reads_mapped"] = bamqc.pct_reads_mapped
stats["bam_med_dp"] = bamqc.med_dp
stats["bam_depth_10"] = bamqc.genome_cov[10]
stats["bam_depth_5"] = bamqc.genome_cov[20]

json.dump(stats,open("%s.stats.json" % prefix,"w"))
O = open("%s.log"%prefix,"w")
for x in ["fastq_mean_read_len","fastq_read_num","kraken_pct_filt","bam_pct_reads_mapped","bam_med_dp","bam_depth_5","bam_depth_10"]:
	O.write("%s\t%s\n" % (x,stats[x]))
O.close()
