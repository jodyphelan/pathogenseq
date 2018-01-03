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
	def __init__(self,bam_file,prefix,ref_file=None):
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
