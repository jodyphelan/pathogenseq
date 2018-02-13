from __future__ import division
import sys
import subprocess
import json
import os
from files import *
from fasta import *
from mvcf import *


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def create_mappability_file(ref_file,threads):
	cmd = "gem-indexer -i %s -o genome" % ref_file
	run_cmd(cmd)
	cmd = "gem-mappability -I genome.gem -o genome -l 51 -T %s" % threads
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
		O.write("%s\t%s\t%s\t%s\n" % (chrom,lines[i][1],lines[i+1][1]-1,lines[i][2]))
	O.close()

class vcf_merge:
	"""
	A class to manage the merging of gVCF files and subsequence filtering

	Args:
		sample_file(str):
		File containing the sample names with 1 per line
		ref_file(str): Reference file
		prefix(str): Prefix for all output files
		mappability_file(str): Name of the mappability file (will generage if NoneType given)
		vcf_dir(str): Directory contianing the gVCF files
		min_dp(int): Minimum depth for calls
		fmiss(float): The maximum fraction of missing data to keep SNP position
		miss_cut(float): Max fraction cutoff for missing calls per sample
		mix_cut(float): Max fraction cutoff for mixed calls per sample
		low_cov(bool): Optimise filtering for low coverage sequencing
		bed_include(str): Only include variants in regions specified in the BED files
		bed_exclude(str): Exclude variants in regions specified in the BED files
		vcf_ext(str): The extension the the VCF files have (after sample name)
		threads(int): Number of threads to use

	Returns:
		vcf_merge: A vcf_merge class object
	"""

	def __init__(self,sample_file,ref_file,prefix,mappability_file=None,vcf_dir=".",min_dp=10,keep_samples=None,fmiss=0.1,miss_cut=0.15,mix_cut=0.15,low_cov=False,bed_include=None,bed_exclude=None,threads=4,vcf_ext="vcf.gz"):
		self.params = {}
		self.samples = []
		self.keep_samples = []
		self.params["sample_file"] = sample_file
		self.params["ref_file"] = ref_file
		self.params["threads"] = threads
		self.params["vcf_ext"] = vcf_ext
		self.params["prefix"] = prefix
		self.params["merged_bcf"] = "%s.raw.bcf" % prefix
		self.params["prefilt_bcf"] = "%s.prefilt.bcf" % prefix
		self.params["uniq_filt_bcf"] = "%s.uniq.bcf" % prefix
		self.params["sample_filt_bcf"] = "%s.sample_filt.bcf" % prefix
		self.params["mix_masked_bcf"] = "%s.mix_masked.bcf" % prefix
		self.params["snp_fasta"] = "%s.snps.fasta" % prefix
		self.params["min_dp"] = min_dp
		self.params["fmiss"] = fmiss
		self.params["vcf_dir"] = vcf_dir
		self.params["lq_sample_file"] = "%s.LQ.samples.txt" % prefix
		self.params["hq_sample_file"] = "%s.HQ.samples.txt" % prefix
		self.params["miss_cut"] = miss_cut
		self.params["mix_cut"] = mix_cut
		self.params["low_cov"] = low_cov
		self.params["qual_file"] = "%s.sample_quals.txt" % prefix
		self.params["bed_include"] = "bcftools view -T %s |" % bed_include if bed_include!=None else ""
		self.params["bed_exclude"] = "bcftools view -T ^%s |" % bed_exclude if bed_exclude!=None else ""
		if not mappability_file:
			create_mappability_file(ref_file,threads)
			self.params["mappability_file"] = "genome.mappability.bed"
		else:
			filecheck(mappability_file)
			self.params["mappability_file"] = mappability_file
		for l in open(sample_file):
			self.samples.append(l.rstrip())
		for s in self.samples:
			self.params["temp"] = s
			filecheck("%(vcf_dir)s/%(temp)s.%(vcf_ext)s" %self.params)
		if keep_samples and filecheck(keep_samples):
			self.keep_samples = [x.rstrip() for x in open(keep_samples).readlines()]
		filecheck(sample_file)
		filecheck(ref_file)


	def merge(self):
		"""Merge gVCF files"""
		cmd = "cat %(sample_file)s | xargs -i -P %(threads)s sh -c \"if [ ! -f %(vcf_dir)s/{}.%(vcf_ext)s.csi ]; then bcftools index %(vcf_dir)s/{}.%(vcf_ext)s; fi;\"" % self.params
		run_cmd(cmd)
		tmp_bcfs = []
		for i,tmp_samples in enumerate(list(chunks(self.samples,500))):
			self.params["tmp_bcf"] = "%s.%s.tmp.bcf" % (prefix,i)
			self.params["vcf_files"] = " ".join(["%s/%s.%s" % (self.params["vcf_dir"],x,self.params["vcf_ext"]) for x in tmp_samples])
			cmd = "bcftools merge --threads %(threads)s -g %(ref_file)s -o %(tmp_bcf)s -O b %(vcf_files)s" % self.params
			run_cmd(cmd)
			tmp_bcfs.append(self.params["tmp_bcf"])
		self.params["vcf_files"] = " ".join(tmp_bcfs)
		cmd = "bcftools merge --threads %(threads)s -g %(ref_file)s -o %(tmp_bcf)s -O b %(vcf_files)s" % self.params
		run_cmd(cmd)


	def extract_variants(self):
		"""Extract all variant positions"""
		cmd = "bcftools +setGT %(merged_bcf)s -- -t q -i 'FMT/DP<%(min_dp)s' -n . | %(bed_include)s %(bed_exclude)s bcftools view -i 'AC>=0 && F_MISSING<%(fmiss)s' -o %(prefilt_bcf)s -O b" % self.params
		run_cmd(cmd)

	def filt_non_uniq(self):
		"""Filter out non unique positions"""
		non_uniq = []
		O = open("genome.non_uniq.bed","w")
		for l in open(self.params["mappability_file"]):
			arr = l.rstrip().split()
			if float(arr[3])<1:
				O.write(l)
		O.close()
		self.params["non_uniq_bed"] = "genome.non_uniq.bed"
		cmd = "bcftools view -T ^%(non_uniq_bed)s %(prefilt_bcf)s -O b -o %(uniq_filt_bcf)s" % self.params
		run_cmd(cmd)

	def sample_filt(self):
		"""Filter out low quality samples"""
		num_calls = int(subprocess.Popen("bcftools view %(uniq_filt_bcf)s -H | wc -l" % self.params,shell=True,stdout=subprocess.PIPE).communicate()[0].rstrip())
		cmd = "cat %(sample_file)s | xargs -i -P %(threads)s sh -c \"bcftools view -s {} %(uniq_filt_bcf)s | bcftools view -g miss -H | wc -l > {}.miss\"" % self.params
		run_cmd(cmd)
		cmd = "cat %(sample_file)s | xargs -i -P %(threads)s sh -c \"bcftools view -s {} %(uniq_filt_bcf)s | bcftools view -g het -H | wc -l > {}.mix\"" % self.params
		run_cmd(cmd)
		miss = {}
		mix = {}
		self.lq_samples = []
		self.hq_samples = []
		HQ = open(self.params["hq_sample_file"],"w")
		LQ = open(self.params["lq_sample_file"],"w")
		QF = open(self.params["qual_file"],"w")
		QF.write("sample\tmix\tmiss\n")
		for s in self.samples:
			miss[s] = int(open("%s.miss" % s).readline().rstrip())/num_calls
			mix[s] = int(open("%s.mix" % s).readline().rstrip())/num_calls
			QF.write("%s\t%s\t%s\n" % (s,mix[s],miss[s]))
			os.remove("%s.miss" % (s))
			os.remove("%s.mix" % (s))
			if s in self.keep_samples:
				self.hq_samples.append(s)
				HQ.write("%s\n" % s)
			elif miss[s]>self.params["miss_cut"] or mix[s]>self.params["mix_cut"]:
				self.lq_samples.append(s)
				LQ.write("%s\n" % s)
			else:
				self.hq_samples.append(s)
				HQ.write("%s\n" % s)
		HQ.close()
		LQ.close()
		QF.close()
		cmd = "bcftools view -S %(hq_sample_file)s -a -c 1 -o %(sample_filt_bcf)s -O b %(uniq_filt_bcf)s" % self.params
		run_cmd(cmd)

	def mask_mixed(self):
		"""Create a BCF file with mixed called masked as missing"""
		cmd = "bcftools +setGT %(sample_filt_bcf)s -- -t q -i 'GT=\"het\"' -n . | bcftools view -Ob -o %(mix_masked_bcf)s" % self.params
		run_cmd(cmd)
	def get_bcf_obj(self):
		return bcf(self.params["mix_masked_bcf"],self.params["ref_file"])
