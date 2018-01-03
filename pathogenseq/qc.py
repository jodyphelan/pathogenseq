from __future__ import division
import sys
import subprocess
from files import *
import numpy as np
import gzip
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import json
import re
from fasta import *
from fastq import *
from bam import *
################################
########## Functions ###########
################################

def gsize_convert(x):
	d = {"G":1e9,"M":1e6,"K":1e3}
	num = float(x[:-1])
	char = x[-1]
	if char not in d: print "%s not a valid value";quit()
	return num*d[char]


def get_genome_cov(bam_file,ref_file,min_dp):
	fdict = fasta(ref_file).fa_dict
	ref_cov = {}
	for s in fdict:
		ref_cov[s] = [0 for x in range(len(fdict[s]))]
	samtools_cmd = "samtools depth -aa --reference %s %s" % (ref_file,bam_file)
	for l in subprocess.Popen(samtools_cmd,shell=True,stdout=subprocess.PIPE).stdout:
		arr = l.rstrip().split()
		if arr[0] not in ref_cov: print "Can't find %s in FASTA...Have you used the correct reference sequence?";quit()
		ref_cov[arr[0]][int(arr[1])-1] = int(arr[2])
	all_genome = []
	for s in fdict:
		all_genome+=ref_cov[s]
	genome_cov = {}
	for dp in min_dp:
		genome_cov[dp] = len([1 for d in all_genome if d>dp])/len(all_genome)

	med = int(np.median(all_genome))
	return genome_cov,med,ref_cov

def flagstat(bam_file):
	lines = []
	samtools_cmd = "samtools flagstat %s" % (bam_file)
	for l in subprocess.Popen(samtools_cmd,shell=True,stdout=subprocess.PIPE).stdout:
		arr = l.rstrip().split()
		lines.append(arr)
	return float(lines[4][4][1:-1])


################################
########### Classes ############
################################

class qc_fastq:
	"""
	A class to extract basic QC stats and run QC programs on fastQ files

	Args:
		prefix(str): Prefix for output files
		fq1(str): First read file [required]
		fq2(str): Second read file. Pass NoneType if there is no second read
		optimise(Bool): Choose if you want to calculate metrics based on whole read file or based on the first 10 percent.
		threads(int): Number of threads to use for multithreaded methods
		kraken_db(str):	Location of the kraken database (if needed)

	Returns:
		qc_fastq: A qc_fastq class object
	"""
	params = {"fq1":"","fq2":""}
	read_len = []
	read_num = 0
	paired = False
	kraken_run = False
	def __init__(self,prefix,fq1,fq2=None,optimise=True,threads=20):
		if filecheck(fq1):
			self.params["fq1"] = fq1
		if fq2 and filecheck(fq2):
			self.params["fq2"] = fq2
			self.paried = True
		self.params["prefix"] = prefix
		self.params["threads"] = threads


		self.read_num = int(gz_file_len(fq1)/4)*2 if self.paired else int(gz_file_len(fq1)/4)
		self.read_pairs = self.read_num/2 if self.paired else self.read_num
		FQ = gzip.open(fq1)
		i = int(self.read_pairs*0.10) if optimise else self.read_pairs
		for j in range(i):
			FQ.readline()
			self.read_len.append(len(FQ.readline()))
			FQ.readline()
			FQ.readline()
		self.median_read_len = np.median(self.read_len)
		self.mean_read_len = np.mean(self.read_len)
	def approx_depth(self,genome_size):
		"""Return approx depth for a given genome size"""
		return self.read_num*self.mean_read_len/gsize_convert(genome_size)
	def run_kraken(self,kraken_db,filter_fastq = None):
		"""
		Run kraken with an option to create filtered fastq files

		Args:
			filter_fastq(str): NCBI Taxonomy code use when extracting reads
		"""
		self.params["kraken_db"] = kraken_db
		self.params["kraken_file"] = "%(prefix)s.kraken" % self.params
		cmd = "kraken --db %(kraken_db)s --threads %(threads)s --fastq-input --gzip-compressed --output %(kraken_file)s --paired --check-names %(fq1)s %(fq2)s" % self.params
		run_cmd(cmd)

		if filter_fastq:
			taxa = filter_fastq.split(",")
			o1 = "%(prefix)s_1.kraken_filt.fastq.gz" % self.params
			o2 = "%(prefix)s_2.kraken_filt.fastq.gz" % self.params
			self.params["kr_filt_fq_1"] = o1
			self.params["kr_filt_fq_2"] = o2
			readnames = set()
			for l in open(self.params["kraken_file"]):
				arr = l.rstrip().split()
				if arr[2] in taxa:
					readnames.add(arr[1])
			R1 = gzip.open(self.params["fq1"])
			R2 = gzip.open(self.params["fq2"])
			O1 = gzip.open(o1,"wb")
			O2 = gzip.open(o2,"wb")

			for seqname1 in R1:
				seqname1 = seqname1.rstrip()
				seq1 = next(R1).rstrip()
				next(R1)
				qual1 = next(R1).rstrip()
				seqname2 = next(R2).rstrip()
				seq2 = next(R2).rstrip()
				next(R2)
				qual2 = next(R2).rstrip()
				if seqname1.split()[0][1:] in readnames:
					O1.write("%s\n%s\n+\n%s\n" % (seqname1,seq1,qual1))
					O2.write("%s\n%s\n+\n%s\n" % (seqname2,seq2,qual2))
			O1.close()
			O2.close()
			self.run_kraken = True
		def get_mapper_from_kraken(self,ref):
			"""
			Get a mapping class from the kraken filtered fastQ files

			Returns:
				mapping: A mapping class object
			"""
			if self.kraken_run==False:
				print "Please run kraken filtering first...exiting"
				quit()
			return mapping(self.params["kr_filt_fq_1"],self.params["kr_filt_fq_2"],ref,self.prefix,threads=20,platform="Illumina",call_method="optimise",mapper="bwa")



class qc_bam:
	"""
	A class to extract basic QC stats and run QC programs on fastQ files

	Args:
		bam(str): Bam file
		ref(str): Refrence fasta
		cov_thresholds(list): List of integers to use in the percentage genome covered calculation

	Returns:
		qc_bam: A qc_bam class object
	"""
	bam = None
	ref = None
	def __init__(self,bam,ref,cov_thresholds=[5,10,20]):
		if filecheck(bam): self.bam = bam
		if filecheck(ref): self.ref = ref
		self.genome_cov,self.med_dp,self.ref_dp = get_genome_cov(bam,ref,cov_thresholds)
		self.pct_reads_mapped = flagstat(bam)
	def plot_cov(self,chrom,imgfile,window=10000,step=5000,optimise=True):
		"""
		Plot coverage across chromosomes

		Args:
			chrom(str):	Chromosome name
			imgfile(str): Name of the output png
			window(int): Window size for the sliding window coverage calculation
			step(int): Step size for the sliding window coverage calculation
			optimise(bool): Optimise window and step size for chromosome len
		"""
		chrom_size = len(self.ref_dp[chrom])
		if chrom_size<100000:
			n,d = "K",1000
		elif chrom_size>100000 and chrom_size<1000000000:
			n,d = "M",1000000
		else:
			n,d = "G",1000000000
		if optimise:
			if chrom_size<100000:
				window,step=100,50
			elif chrom_size>100000 and chrom_size<1000000:
				window,step=1000,500

		x = []
		y = []
		hw = int(window/2)
		for i in range(hw,len(self.ref_dp[chrom])-hw,step):
			x.append(i/d)
			y.append(int(np.median(self.ref_dp[chrom][i-hw:i+hw])))
		fig = plt.figure()
		plot = fig.add_subplot(111)
		plot.plot(x,y)
		plot.set_ylim(bottom=0)
		if max(y)>200:
			plot.set_yscale('symlog')
		plot.set_xlabel("Genome Position (%sb)" % n)
		plot.set_ylabel("Median Coverage (Window size:%s)" % window)
		fig.savefig(imgfile)
	def save_cov(self,filename):
		"""Save coverage to a json file"""
		json.dump(self.ref_dp,open(filename,"w"))
	def region_cov(self,regions):
		"""
		Return a dictionary with mean depth across selected regions

		Args:
			regions(list): A list with each element consisting of a ``tuple`` with 4 strings: 1) chromosome, 2) start, 3) end and 4) ID

		Returns:
			dict: A dictionary with mean depth across selected regions
		"""
		results = {}
		for chrom,start,end,name in regions:
			results[name] = np.mean(self.ref_dp[chrom][int(start)-1:int(end)])
		return results
	def bed_cov(self,bed_file):
		"""
		Return a dictionary with mean depth across selected regions in BED file

		Args:
			bed_file(str): A bed file with the 4th column containing the region ID

		Returns:
			dict: A dictionary with mean depth across selected regions
		"""

		regions = []
		for l in open(bed_file):
			#Chromosome	start	end	name
			arr = l.rstrip().split()
			regions.append(tuple(arr[:4]))
		return self.region_cov(regions)
	def gff_cov(self,gff_file,key="ID"):
		"""
		Return a dictionary with mean depth across selected regions in GFF file

		Args:
			gff_file(str): A gff file region coordinates
			key(str): A key to use as the region ID (e.g. for ID=katG the key is 'ID')

		Returns:
			dict: A dictionary with mean depth across selected regions
		"""
		regions = []
		key_re = re.compile("%s=([\w\.\-\_]+)"%key)
		for l in open(gff_file):
			if l[0]=="#": continue
			arr = l.rstrip().split()
			if "%s="%key not in l:
				print "Warining: %s not found in %s" % (key,l)
				continue
			name = key_re.search(l).group(1)
			regions.append((arr[0],arr[3],arr[4],name))
		return self.region_cov(regions)
