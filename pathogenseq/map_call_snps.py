from __future__ import division
import subprocess
from files import *
import os
import numpy as np
from qc import *

class fastq:
	"""
	Class for performing mapping to reference.

	Args:
		r1(str): Location of the first read [required]
		r2(str): Location of the second read (if any) use NoneType if there is none.
		ref_file(str): Location of the refrence file [required]
		prefix(str): Prefix for the output files generated (e.g. <prefix>.bam)
		threads(int): Number of threads to use
		platform(str): Platform used in sequencing (currently only 'Illumina' is supported)
		call_method(str): Method used to call variants. Can be ``optimise``, ``low`` or ``high``.
		mapper(str): Mapping tool. Can be ``bwa`` or ``bowtie``

	Returns:
		mapping: A mapping class object
	"""
	params = {}
	paired = False
	call_method = "high"
	mapper = "bwa"
	def __init__(self,r1,r2,ref_file,prefix,threads=20,platform="Illumina",mapper="bwa"):
		if r1 and filecheck(r1):
			self.params["r1"] = r1
		else:
			print "Provide at least one fastq file...Exiting\n";quit()
		if r2 and filecheck(r2):
			self.params["r2"] = r2
			self.paired = True
		self.mapper = mapper
		self.params["prefix"] = prefix
		self.params["threads"] = threads
		self.params["platform"] = platform
		if filecheck(ref_file):
			self.params["ref_file"] = ref_file
		if nofile(ref_file+".bwt"):
			bwa_index(ref_file)
		self.params["bam_file"] = "%s.bam" % prefix
		self.params["r1_tp"] = "%s_1_trim_pair.fq" % prefix
		self.params["r1_tu"] = "%s_1_trim_unpair.fq" % prefix
		self.params["r2_tp"] = "%s_2_trim_pair.fq" % prefix
		self.params["r2_tu"] = "%s_2_trim_unpair.fq" % prefix
		self.params["r1t"] = "%s_trim.fq" % prefix
		self.call_method = call_method

	def trim(self):
		"""Perform trimming"""
		if self.paired:
			cmd = "java -jar ~/bin/trimmomatic.jar PE -threads %(threads)s -phred33 %(r1)s %(r2)s %(r1_tp)s %(r1_tu)s %(r2_tp)s %(r2_tu)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % self.params
		else:
			cmd = "java -jar ~/bin/trimmomatic.jar SE -threads %(threads)s -phred33 %(r1)s %(r1t)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % self.params
		run_cmd(cmd)

	def map(self):
		"""Perform mapping"""
		prefix = self.params["prefix"]
		self.params["bwa_prefix"] = "bwa mem -t %(threads)s -c 100 -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' -M -T 50" % self.params
		if self.paired:
			self.params["bam_pair"] = "%s.pair.bam" % prefix
			self.params["bam_single1"] = "%s.single1.bam" % prefix
			self.params["bam_single2"] = "%s.single2.bam" % prefix
			self.params["bam_unsort"] = "%s.unsort.bam" % prefix
			self.params["temp"] = "%s.paired.bam" % self.params["prefix"]
			if self.mapper=="bwa":
				cmd = "%(bwa_prefix)s %(ref_file)s %(r1_tp)s %(r2_tp)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair)s -" % self.params
				run_cmd(cmd)
				cmd = "%(bwa_prefix)s %(ref_file)s %(r1_tu)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single1)s" % self.params
				run_cmd(cmd)
				cmd = "%(bwa_prefix)s %(ref_file)s %(r2_tu)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single2)s" % self.params
				run_cmd(cmd)
			elif self.mapper=="bowtie":
				cmd = "bowtie2 -p %(threads)s -x %(ref_file)s -1 %(r1_tp) -2 %(r2_tp) | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair)s -" % self.params
				run_cmd(cmd)
				cmd = "bowtie2 -p %(threads)s -x %(ref_file)s -U %(r1_tu)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single1)s" % self.params
				run_cmd(cmd)
				cmd = "bowtie2 -p %(threads)s -x %(ref_file)s -U %(r2_tu)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single1)s" % self.params
				run_cmd(cmd)
			cmd = "samtools merge -@ %(threads)s %(bam_unsort)s %(bam_pair)s %(bam_single1)s %(bam_single2)s" % self.params
			run_cmd(cmd)
			cmd = "samtools sort -@ %(threads)s -o %(bam_file)s %(bam_unsort)s" % self.params
			run_cmd(cmd)
			rm_files([self.params["r1_tp"],self.params["r1_tu"],self.params["r2_tp"],self.params["r2_tu"],self.params["bam_pair"],self.params["bam_single1"],self.params["bam_single2"],self.params["bam_unsort"]])
		else:
			if self.mapper=="bwa":
				cmd = "%(bwa_prefix)s %(ref_file)s %(r1t)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % self.params
			elif self.mapper=="bowtie":
				cmd = "bowtie2 -p %(threads)s -x %(ref_file)s -U %(r1t)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s" % self.params
			run_cmd(cmd)
			rm_files([self.params["r1t"]])



class bam:
	params = {}
	def __init__(self,bam_file,prefix,ref_file=None):
		if filecheck(bam_file): self.params["bam_file"] = bam_file
		self.params["prefix"] = prefix
		if ref_file:
			if filecheck(ref_file): self.params["ref_file"] = ref_file
		else:
			self.params["ref_file"] = ref_file

	def call_snps(self,ref_file=None,call_method="optimise"):
		"""
		Perform SNP calling

		Args:
			ref_file(str): reference file (not required if passed to the bam initiator).
			call_method(str): optimise variant calling based on high or low depth. Options: high|low|optimise
		"""
		self.params["vcf_file"] = "%s.vcf.gz" % self.params["prefix"]
		if ref_file:
			if filecheck(ref_file): self.params["ref_file"] = ref_file
		else:
			if self.params["ref_file"]==None and ref_file==None:
				print "Please provide a reference fasta file...Exiting"
				quit()
		if call_method=="optimise":
			dp = []
			cmd = "samtools depth %(bam_file)s" % self.params
			print "Optimising call method"
			for l in subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout:
				arr = l.rstrip().split()
				dp.append(int(arr[2]))
			med_dp = np.median(dp)
			print "Median depth: %s" % med_dp
			if med_dp<30:
				print "Using low depth approach"
				call_method = "low"
			else:
				print "Using high depth approach"
				call_method = "high"

		if call_method=="high":
			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -t DP | bcftools call -mg 10 -V indels -Oz -o %(vcf_file)s" % self.params
		else:
			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -ABq0 -Q0 -t DP | bcftools call -mg 10 -V indels -Oz -o %(vcf_file)s" % self.params
		run_cmd(cmd)
	def get_bam_qc(self,ref_file,cov_thresholds=[1,5,10,20]):
		"""
		Get a qc_bam object

		Args:
			cov_thresholds(list): List of integers to use in the percentage genome covered calculation
		Returns:
			qc_bam: A qc_bam object
		"""
		if ref_file:
			if filecheck(ref_file): self.params["ref_file"] = ref_file
		else:
			if ref_file==None:
				print "Please provide a reference fasta file...Exiting"
				quit()
		return qc_bam(self.params["bam_file"],self.params["ref_file"],cov_thresholds)
