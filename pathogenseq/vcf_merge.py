from __future__ import division
import sys
import subprocess
import json
import os
from .files import *
from .fasta import *
from .mvcf import *


def split_list(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def create_mappability_file(ref_file,threads):
	cmd = "gem-indexer -i %s -o genome" % ref_file
	run_cmd(cmd)
	cmd = "gem-mappability -I genome.gem -o genome -l 65 -T %s" % threads
	run_cmd(cmd)
	cmd = "gem-2-wig -I genome.gem -i genome.mappability -o genome"
	run_cmd(cmd)

	lines = []
	for l in open("genome.wig"):
		arr = l.rstrip().split()
		if arr[0]=="variableStep":
			chrom = arr[1].split()[0][6:]
			continue
		lines.append((chrom,int(arr[0]),arr[1]))
	O = open("genome.mappability.bed","w")
	for i in range(len(lines)-1):
		O.write("%s\t%s\t%s\t%s\n" % (lines[i][0],lines[i][1],lines[i+1][1]-1,lines[i][2]))
	O.close()

class vcf_merge:
	"""
	A class to manage the merging of gVCF files and subsequence filtering

	Args:
		sample_file(str):
		File containing the sample names with 1 per line
		ref_file(str): Reference file
		prefix(str): Prefix for all output files
		vcf_dir(str): Directory contianing the gVCF files
		vcf_ext(str): The extension the the VCF files have (after sample name)
		threads(int): Number of threads to use

	Returns:
		vcf_merge: A vcf_merge class object
	"""

	def __init__(self,sample_file,ref_file,prefix,vcf_dir=".",vcf_ext="gbcf",threads=4,min_dp=10):
		add_arguments_to_self(self,locals())
		self.samples = []
		filecheck(sample_file)
		filecheck(ref_file)
		for l in open(sample_file):
			self.samples.append(l.rstrip())
		for s in self.samples:
			self.temp = s
			filecheck("%(vcf_dir)s/%(temp)s.%(vcf_ext)s" % vars(self))
		self.merged_bcf = "%s.raw.bcf" % self.prefix
	def merge(self):
		"""Merge gVCF files"""
		cmd = "cat %(sample_file)s | parallel -j %(threads)s \"if [ ! -f %(vcf_dir)s/{}.%(vcf_ext)s.csi ]; then bcftools index %(vcf_dir)s/{}.%(vcf_ext)s; fi;\"" % vars(self)
		run_cmd(cmd)
		tmp_bcfs = []
		self.tmp_file = "%s.cmd.xargs.txt" % self.prefix
		X = open(self.tmp_file,"w")
		chunk_size = 100
		chunks = list(split_list(self.samples,chunk_size))
		if len(chunks[-1])==1:
			chunks[0] = chunks[0] + chunks.pop()
		for i,tmp_samples in enumerate(chunks):
			self.tmp_bcf = "%s.%s.tmp.bcf" % (self.prefix,i)
			self.vcf_files = " ".join(["%s/%s.%s" % (self.vcf_dir,x,self.vcf_ext) for x in tmp_samples])
			cmd = "bcftools merge --threads 2 -O u -g %(ref_file)s %(vcf_files)s | bcftools filter -e 'FMT/DP<%(min_dp)s' -o %(tmp_bcf)s -O b -S . && bcftools index --threads 2 %(tmp_bcf)s" % vars(self)
			X.write("%s\n"%cmd)
			tmp_bcfs.append(self.tmp_bcf)
		X.close()
		tmp_threads = 1 if self.threads==1 else int(self.threads/2)
		if tmp_threads>8: tmp_threads = 8 #Stop the max number of files being hit... assuming 1000 files is limit
		cmd = "cat %s | parallel -j %s" % (self.tmp_file,tmp_threads)
		run_cmd(cmd)
		if len(tmp_bcfs)>1:
			self.vcf_files = " ".join(tmp_bcfs)
			cmd = "bcftools merge --threads %(threads)s -g %(ref_file)s -o %(merged_bcf)s -O b %(vcf_files)s" % vars(self)
			run_cmd(cmd)
			rm_files(tmp_bcfs)
		else:
			cmd = "mv %(tmp_bcf)s %(merged_bcf)s" % vars(self)
			run_cmd(cmd)
		return bcf(self.merged_bcf)
