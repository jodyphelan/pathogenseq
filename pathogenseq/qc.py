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

################################
########## Functions ###########
################################

def gsize_convert(x):
	d = {"G":1e9,"M":1e6,"K":1e3}
	num = float(x[:-1])
	char = x[-1]
	if char not in d: print "%s not a valid value";quit()
	return num*d[char]

def fa2dict(filename):
        fa_dict = {}
        seq_name = ""
        for l in open(filename):
                line = l.rstrip()
                if line[0] == ">":
                        seq_name = line[1:].split()[0]
                        fa_dict[seq_name] = []
                else:
                        fa_dict[seq_name].append(line)
        result = {}
        for seq in fa_dict:
                result[seq] = "".join(fa_dict[seq])
	return result

def get_genome_cov(bam_file,ref_file,min_dp):
	fdict = fa2dict(ref_file)
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
	params = {"fq1":"","fq2":""}
	read_len = [] 
	read_num = 0
	paired = False
	def __init__(self,prefix,fq1,fq2=False,method="optimise",threads=20,kraken_db = False):
		if filecheck(fq1):
			self.params["fq1"] = fq1
		if fq2 and filecheck(fq2): 	
			self.params["fq2"] = fq2
			self.paried = True
		self.params["prefix"] = prefix
		self.params["threads"] = threads
		self.params["kraken_db"] = kraken_db

		self.read_num = int(gz_file_len(fq1)/4)*2 if self.paired else int(gz_file_len(fq1)/4)
		self.read_pairs = self.read_num/2 if self.paired else self.read_num
		FQ = gzip.open(fq1)
		i = int(self.read_pairs*0.10) if method=="optimise" else self.read_pairs
		for j in range(i):
			FQ.readline()
			self.read_len.append(len(FQ.readline()))
			FQ.readline()
			FQ.readline()
		self.median_read_len = np.median(self.read_len)
		self.mean_read_len = np.mean(self.read_len)
	def approx_depth(self,genome_size):
		return self.read_num*self.mean_read_len/gsize_convert(genome_size)
	def run_kraken(self,filter_fastq = False):
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
				


	
class qc_bam:
	bam = None
	ref = None
	def __init__(self,bam,ref,cov_thresholds=[5,10,20]):
		if filecheck(bam): self.bam = bam
		if filecheck(ref): self.ref = ref
		self.genome_cov,self.med_dp,self.ref_dp = get_genome_cov(bam,ref,cov_thresholds)
		self.pct_reads_mapped = flagstat(bam)
	def plot_cov(self,chrom,imgfile,window=10000,step=5000,method="optimise"):
		chrom_size = len(self.ref_dp[chrom])
		if chrom_size<100000:
			n,d = "K",1000
		elif chrom_size>100000 and chrom_size<1000000000:
			n,d = "M",1000000
		else:
			n,d = "G",1000000000
		if method=="optimise":
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
		json.dump(self.ref_dp,open(filename,"w"))
	def region_cov(self,regions):
		results = {}
		for chrom,start,end,name in regions:
			results[name] = np.mean(self.ref_dp[chrom][int(start)-1:int(end)])
		return results
	def bed_cov(self,bed_file):
		regions = []
		for l in open(bed_file):
			#Chromosome	start	end	name
			arr = l.rstrip().split()
			regions.append(tuple(arr))
		return self.region_cov(regions)

