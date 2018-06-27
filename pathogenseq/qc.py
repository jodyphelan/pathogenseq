from __future__ import division
import sys
import subprocess
from .files import *
from .fasta import *
from .fastq import *
import numpy as np
import gzip
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import json
import re

from collections import defaultdict
################################
########## Functions ###########
################################

def gsize_convert(x):
	d = {"G":1e9,"M":1e6,"K":1e3}
	num = float(x[:-1])
	char = x[-1]
	if char not in d: log("%s not a valid value");quit()
	return num*d[char]


def get_genome_cov(bam_file,ref_file,min_dp):
	fdict = fasta(ref_file).fa_dict
	ref_cov = {}
	for s in fdict:
		ref_cov[s] = [0 for x in range(len(fdict[s]))]
	samtools_cmd = "samtools depth -aa --reference %s %s" % (ref_file,bam_file)
	log(samtools_cmd)
	for l in subprocess.Popen(samtools_cmd,shell=True,stdout=subprocess.PIPE).stdout:
		arr = l.rstrip().split()
		if arr[0] not in ref_cov: log("Can't find %s in FASTA...Have you used the correct reference sequence?");quit()
		ref_cov[arr[0]][int(arr[1])-1] = int(arr[2])
	all_genome = []
	for s in fdict:
		all_genome+=ref_cov[s]
	genome_cov = {}
	for dp in min_dp:
		genome_cov[dp] = len([1 for d in all_genome if d>=dp])/len(all_genome)

	med = int(np.median(all_genome))
	return genome_cov,med,ref_cov

def flagstat(bam_file):
	lines = []
	samtools_cmd = "samtools flagstat %s" % (bam_file)
	for l in subprocess.Popen(samtools_cmd,shell=True,stdout=subprocess.PIPE).stdout:
		arr = l.rstrip().split()
		lines.append(arr)
	num = int(lines[4][0])
	pct = 0.0 if num==0 else float(lines[4][4][1:-1])
	return num,pct


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

	def __init__(self,prefix,fq1,fq2=None,optimise=True,threads=4):
		self.params = {"fq1":"","fq2":""}
		self.read_len = []
		self.read_num = 0
		self.paired = False
		self.kraken_run = False
		if filecheck(fq1):
			self.params["fq1"] = fq1
		if fq2 and filecheck(fq2):
			self.params["fq2"] = fq2
			self.paired = True

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
	def run_centrifuge(self,centrifuge_db,filter_fastq=None,threads=4):
		self.params["centrifuge_db"] = centrifuge_db
		self.params["centrifuge_report"] = "%(prefix)s.centrifuge.report.txt" % self.params
		self.params["centrifuge_log"] = "%(prefix)s.centrifuge.log" % self.params
		self.params["threads"] = threads
		if self.paired:
			cmd = "centrifuge -x %(centrifuge_db)s -1 %(fq1)s -2 %(fq2)s -S %(centrifuge_log)s --report-file %(centrifuge_report)s -p %(threads)s" % self.params
		else:
			cmd = "centrifuge -x %(centrifuge_db)s -U %(fq1)s -S %(centrifuge_log)s --report-file %(centrifuge_report)s -p %(threads)s" % self.params
		run_cmd(cmd)
		if filter_fastq:
			num_mtb = 0
			taxa = filter_fastq.split(",")
			self.params["cf_filt_fq_1"] = "%(prefix)s_1.centrifuge_filt.fastq.gz" % self.params
			self.params["cf_filt_fq_2"] = "%(prefix)s_2.centrifuge_filt.fastq.gz" % self.params
			self.params["tmp_file"] = get_random_file()
			read_names = set()
			for l in open(self.params["centrifuge_log"]):
				#K00250:202:HNN53BBXX:8:1101:6066:998    NC_016947.1     1138382 21377   21377   235     302     4
				row = l.rstrip().split()
				if row[2] in taxa:
					num_mtb+=1
					read_names.add(row[0])

			O = open(self.params["tmp_file"],"w")
			O.write("\n".join(list(read_names)))
			O.close()
			cmd = "seqtk subseq %(fq1)s %(tmp_file)s | pigz -p %(threads)s -c > %(cf_filt_fq_1)s" % self.params
			run_cmd(cmd)
			cmd = "seqtk subseq %(fq2)s %(tmp_file)s | pigz -p %(threads)s -c > %(cf_filt_fq_2)s" % self.params
			run_cmd(cmd)
			rm_files([self.params["tmp_file"]])

		top_hit = ""
		top_num_reads = 0
		for l in open(self.params["centrifuge_report"]):
			#Mycobacterium avium     1764    species 6256976 13835   352     2.10431e-06
			row = l.rstrip().split("\t")
			if row[0]=="name": continue
			if int(row[5])>top_num_reads:
				top_hit = row[0].replace(" ","_")
				top_num_reads = int(row[4])

		if filter_fastq:
			tmp = [top_hit,num_mtb/self.read_num]
			return self.params["cf_filt_fq_1"],self.params["cf_filt_fq_2"],tmp
		else:
			return top_hit,top_num_reads

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

	def __init__(self,bam,ref,cov_thresholds=[1,5,10],threads=4):
		self.bam = None
		self.ref = None
		self.threads = threads
		if filecheck(bam): self.bam = bam
		if filecheck(ref): self.ref = ref
		self.genome_cov,self.med_dp,self.ref_dp = get_genome_cov(bam,ref,cov_thresholds)
		self.num_reads_mapped,self.pct_reads_mapped = flagstat(bam)

	def plot_cov(self,chrom,imgfile,start=None,end=None,window=10000,step=5000,optimise=True,plot_median=True,primers=None):
		"""
		Plot coverage across chromosomes

		Args:
			chrom(str):	Chromosome name
			imgfile(str): Name of the output png
			window(int): Window size for the sliding window coverage calculation
			step(int): Step size for the sliding window coverage calculation
			optimise(bool): Optimise window and step size for chromosome len
		"""
		if plot_median:
			chrom_med_dp = np.median(self.ref_dp[chrom])
		if start and end:
			region_size = end-start
			offset = int(region_size*0.05)
			new_start = start-offset
			new_end = end+offset
		else:
			offset=False
			region_size = len(self.ref_dp[chrom])
			start = 0
			end = region_size
			new_start = start
			new_end = end
		if region_size<100000:
			n,d = "K",1000
		elif region_size>100000 and region_size<1000000000:
			n,d = "M",1000000
		else:
			n,d = "G",1000000000
		if optimise:
			if region_size<10000:
				window,step=2,1
			elif region_size<100000:
				window,step=100,50
			elif region_size>100000 and region_size<1000000:
				window,step=1000,500
		else:
			if region_size<10000:
				window,step=2,1
		log("Outputting coverage plot for region (%sbp) with window=%s and step=%s" % (region_size,window,step))
		x = []
		y = []
		hw = int(window/2)
		for i in range(new_start+hw,new_end-hw,step):
			x.append(i/d)
			y.append(int(np.median(self.ref_dp[chrom][i-hw:i+hw+1])))
		fig = plt.figure()
		plot = fig.add_subplot(111)
		plot.plot(x,y)
		plot.set_ylim(bottom=0)
		if max(y)>200:
			plot.set_yscale('symlog')
		plot.set_xlabel("Genome Position (%sb)" % n)
		plot.set_ylabel("Median Coverage (Window size:%s)" % window)
		if plot_median:
			ymax = max(y) if max(y)>chrom_med_dp else chrom_med_dp
			plot.set_ylim(top=ymax+ymax*0.05)
			plot.axhline(xmin=0,xmax=1,y=chrom_med_dp,color="orange",linestyle="dashed")
		if offset:
			plot.axvline(ymin=0,ymax=0.05,x=start/d,color="orange")
			plot.axvline(ymin=0,ymax=0.05,x=end/d,color="orange")
		if primers:
			locations = fasta(self.ref).find_primer_positions(primers)
			for primer in sorted(locations,key=lambda x:locations[x]["start"]):
				p = locations[primer]
				plot.plot((p["start"]/d,p["end"]/d),(0,0),'r-',lw=3)
		fig.savefig(imgfile)
	def save_cov(self,filename,bed=None):
		"""Save coverage to a json file"""
		if bed:
			bed_regions = load_bed(bed,[1,2,3],4)
			bed_cov = {d:[] for d in bed_regions.keys()}
			for locus in bed_regions:
				tmp = bed_regions[locus]
				for i in range(int(tmp[1]),int(tmp[2])):
					bed_cov[locus].append(self.ref_dp[tmp[0]][i])
			json.dump(bed_cov,open(filename,"w"))
		else:
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
				log("Warining: %s not found in %s" % (key,l))
				continue
			name = key_re.search(l).group(1)
			regions.append((arr[0],arr[3],arr[4],name))
		return self.region_cov(regions)
	def extract_gc_skew(self,filename,window=1000,step=500):
		fa_dict = fasta(self.ref).fa_dict
#		gc = []
#		cov = []
		hw = int(window/2)
		results = defaultdict(list)
		for s in fa_dict:
			for i in range(hw,len(fa_dict[s])-hw,step):
				seq = fa_dict[s][i-hw:i+hw]
				tmp = dict((c, seq.count(c)) for c in ["C","G"])
				results[int((tmp["G"]+tmp["C"])/(window)*100)].append(int(np.median(self.ref_dp[s][i-hw:i+hw])))
#				gc.append(int((tmp["G"]+tmp["C"])/(window)*100))
#				cov.append(int(np.median(self.ref_dp[s][i-hw:i+hw])))
#		O = open(filename,"w")
#		for i in range(len(gc)):
#			O.write("%s\t%s\n" % (gc[i],cov[i]))
#		O.close()
		json.dump(results,open(filename,"w"))
