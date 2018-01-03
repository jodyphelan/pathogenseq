from files import *
from collections import defaultdict


def recode_indels(indels):
	#["C+5CGGGG","G-1C"]
	sorted_indels = sorted([x for x in indels],key=lambda y:len(y))
	largest_del = sorted([x for x in indels if "-" in x],key=lambda y:len(y))

	if len(largest_del)==0:
		leftseq = indelre.search(sorted_indels[-1]).group(1)
		rightseq = ""
	else:
		leftseq = indelre.search(largest_del[-1]).group(1)
		rightseq = indelre.search(largest_del[-1]).group(3)
	refseq = leftseq+rightseq
	recoded_indels = []
	for i in indels:
		if "-" in i:
			indel_len = int(indelre.search(i).group(2))
			iseq = leftseq+rightseq[:-(len(rightseq)-indel_len)]
		elif "+" in i:
			iseq = leftseq+indelre.search(i).group(3)+rightseq
		else:
			iseq = i+rightseq
		recoded_indels.append(iseq)
	return (refseq,recoded_indels)





class bam:
	"""
	A class to perform operations on BAM files such as SNP calling

	Args:
		bam_file(str): The BAM file [required]
		prefix(str): A prefix for output files [required]
		ref_file(ref_file): A reference (needed by some methods)

	Returns:
		bam: A bam class object
	"""
	params = {}
	def __init__(self,bam_file,prefix,ref_file=None,platform="Illumina"):
		if filecheck(bam_file): self.params["bam_file"] = bam_file
		self.params["prefix"] = prefix
		if ref_file:
			if filecheck(ref_file): self.params["ref_file"] = ref_file
		else:
			self.params["ref_file"] = ref_file

	def call_snps(self,ref_file=None,call_method="optimise",min_dp=10):
		"""
		Create a gVCF file (for a description see:https://sites.google.com/site/gvcftools/home/about-gvcf)

		Args:
			ref_file(str): reference file (not required if passed to the bam initiator).
			call_method(str): optimise variant calling based on high or low depth. Options: high|low|optimise
			min_dp(int): Minimum depth required to group site into reference-block
		"""
		self.params["min_dp"] = min_dp
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
			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -t DP | bcftools call -mg %(min_dp)s -V indels -Oz -o %(vcf_file)s" % self.params
		else:
			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -ABq0 -Q0 -t DP | bcftools call -mg %(min_dp)s -V indels -Oz -o %(vcf_file)s" % self.params
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

	def do_pileup(self,bed_file=None):
		self.params["temp"] = bed_file
		self.params["temp_pileup"] = "%(prefix)s.temp.pileup" % self.params
		if bed_file:
			cmd = "samtools view -@ %(threads)s -bL %(temp)s %(bam_file)s  > %(temp_bam)s && htsbox pileup -f %(ref_file)s -Q 8 %(temp_bam)s > %(temp_pileup)s" % self.params
		else:
			cmd = "htsbox pileup -f %(ref_file)s -Q 8 %(bam_file)s > %(temp_pileup)s" % self.params
		run_cmd(cmd)
	def htsbox_calls(self,bed_file):
		bed_pos = set()
		for l in open(bed_file):
			arr = l.rstrip().split()
			for i in range(int(arr[1]),int(arr[2])+1):
				bed_pos.add((arr[0],str(i)))
		self.do_pileup(bed_file=bed_file)
		final_calls = defaultdict(lambda :defaultdict(list))
		if self.params["platform"] == "Illumina":
			for l in open(self.params["temp_pileup"]):
				arr = l.rstrip().split()
				if (arr[0],arr[1]) not in bed_pos: continue
				calls = arr[3].split(",")
				cov = [int(x) for x in arr[4].split(":")[1].split(",")]
				tot = sum(cov)
				for i in range(len(calls)):
					final_calls[arr[0]][arr[1]].append((calls[i],cov[i]/tot,cov[i]))
		elif self.params["platform"] == "minION":
			for l in open(self.params["temp_pileup"]):
				#Chromosome	  23	  G	   G,G-1C,G+3AAA   0/1:49,1,1
				arr = l.rstrip().split()
				if (arr[0],arr[1]) not in bed_pos: continue
				alleles = arr[3].split(",")
				depth = [int(x) for x in arr[4].split(":")[1].split(",")]
				max_allele_dp = max(depth)
				max_allele = alleles[depth.index(max_allele_dp)]
				max_allele_frac = max_allele_dp/sum(depth)
				if len(max_allele)>1:
					max_allele = recode_indels([arr[1],max_allele])[1][0]
				if sum(depth)<min_dp:
					call = "N"
				if max_allele_frac<min_frac:
					call = "N"
				else:
					call = max_allele
				final_calls[arr[0]][arr[1]].append((call,1,max_allele_dp))
		return final_calls
