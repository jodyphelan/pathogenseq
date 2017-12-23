from __future__ import division
import sys
import pathogenseq
import argparse
import json


stats = {}
fastqqc = ps.qc_fastq(prefix,r1,r2,kraken_db="/opt/storage2/ernest/kraken/kraken/standard_db")
stats["fastq_mean_read_len"] = fastqqc.mean_read_len
stats["fastq_read_num"] = fastqqc.read_num

fastqqc.run_kraken("77643,1773,78331,33894,1765")

fr1 = "%s_1.kraken_filt.fastq.gz" %prefix
fr2 = "%s_2.kraken_filt.fastq.gz" %prefix

newfastqqc = ps.qc_fastq(prefix,fr1,fr2)
stats["kraken_pct_filt"] = newfastqqc.read_num/fastqqc.read_num*100

mapper = pathogenseq.mapping(fr1,fr2,ref,prefix,threads=threads)
mapper.trim()
mapper.map()

bamqc = mapper.get_bam_qc()
cov_plot = "%s.cov.png" % (prefix)
bamqc.plot_cov("Chromosome",cov_plot)
stats["bam_pct_reads_mapped"] = bamqc.pct_reads_mapped
stats["bam_med_dp"] = bamqc.med_dp

O = open("%s.log"%prefix,"w")
for x in ["fastq_mean_read_len","fastq_read_num","kraken_pct_filt","bam_pct_reads_mapped","bam_med_dp"]:
	print "%s\t%s" % (x,stats[x])
