from __future__ import division
from .files import *
from .utils import *
from .mvcf import *
from .qc import *
from collections import defaultdict
import re
import numpy as np


import vcf
import pysam
from tqdm import tqdm
def get_overlapping_reads(infile,chrom,start,end,outfile,flank=30,threads=4):
	IN = pysam.AlignmentFile(infile,"rb")
	OUT = pysam.AlignmentFile(outfile,"wb",template=IN)
	if start-flank<0:
		OUT.close()
		return 0
	else:
		i = 0
		for read in IN.fetch(chrom,start,end):
			if read.reference_start<=start-flank and read.reference_end>=end+flank:
				i+=1
				OUT.write(read)
		OUT.close()
		return i



class bam:
	"""
	A class to perform operations on BAM files such as SNP calling

	Args:
		bam_file(str): The BAM file [required]
		prefix(str): A prefix for output files [required]
		ref_file(ref_file): A reference (needed by some methods)
		platform(str): Can be either ``Illumina`` or ``minION``
	Returns:
		bam: A bam class object
	"""
	def __init__(self,bam_file,prefix,ref_file,platform="Illumina",threads=4):
		self.params = {}
		self.bam_file = bam_file
		self.ref_file = ref_file
		index_bam(bam_file,threads=threads)
		if filecheck(bam_file):
			self.params["bam_file"] = bam_file
			self.bam = bam_file
		self.params["prefix"] = prefix
		self.prefix = prefix
		if filecheck(ref_file):
			self.params["ref_file"] = ref_file
			self.ref_fa = fasta(self.params["ref_file"])
			self.ref_fa_dict = self.ref_fa.fa_dict
		self.params["platform"] = platform
		self.params["threads"] = threads
	def generate_primer_bcf(self,threads=4,flank=30):
		self.params["failed_primers"] = "%(prefix)s.failed_primers.bed" % self.params
		primer_ids = []
		FAILED = open(self.params["failed_primers"],"w")
		for l in tqdm(open(self.params["primer_bed_file"])):
			chrom,start,end,pid = l.rstrip().split()[:4]
			primer_ids.append(pid)
			start = int(start)
			end = int(end)
			tmp_bcf = "%s.%s.bcf" % (self.prefix,pid)
			tmp_bam = "%s.%s.bam" % (self.prefix,pid)
			self.params["tmp"] = "%s:%s-%s" % (chrom,start,end)
			self.params["pid"] = pid
			log("Extracting reads for %s" % pid)#
			read_num = get_overlapping_reads(self.bam,chrom,start,end,tmp_bam,flank=30,threads=threads)
			if read_num==0:
				#cmd = "bcftools mpileup  -f %(ref_file)s %(bam_file)s %(mpileup_options)s -r %(tmp)s | bcftools call %(vtype)s -m | bcftools +setGT -Ob -o %(prefix)s.%(pid)s.bcf -- -t a -n ." % self.params
				log("No reads for %s" % pid)
				FAILED.write(l)
			else:
				pass
				#cmd = "samtools index %(prefix)s.%(pid)s.bam && bcftools mpileup  -f %(ref_file)s %(prefix)s.%(pid)s.bam %(mpileup_options)s -r %(tmp)s | bcftools call -t %(tmp)s %(vtype)s -mg %(min_dp)s | bcftools norm -f %(ref_file)s  | bcftools +setGT -Ob -o %(prefix)s.%(pid)s.bcf -- -t q -i 'FMT/DP<%(min_dp)s' -n ." % self.params
		FAILED.close()
		cmd = "cat %(primer_bed_file)s | parallel --progress --col-sep '\\t' -j %(threads)s \"samtools index %(prefix)s.{4}.bam && bcftools mpileup  -f %(ref_file)s %(prefix)s.{4}.bam %(mpileup_options)s -B -r {1}:{2}-{3} | bcftools call -t {1}:{2}-{3} %(vtype)s -mg %(min_dp)s | bcftools norm -f %(ref_file)s  | bcftools +setGT -Ob -o %(prefix)s.{4}.bcf -- -t q -i 'FMT/DP<%(min_dp)s' -n . && bcftools index %(prefix)s.{4}.bcf\"" % self.params
		run_cmd(cmd)
		cmd = "cat %(failed_primers)s | parallel --progress --col-sep '\\t' -j %(threads)s \"bcftools mpileup  -f %(ref_file)s %(bam_file)s %(mpileup_options)s -r {1}:{2}-{3} | bcftools call %(vtype)s -m | bcftools +setGT -Ob -o %(prefix)s.{4}.bcf -- -t a -n .\"" % self.params
		run_cmd(cmd)
		cmd = "bcftools concat `cut -f4 %(primer_bed_file)s | awk '{print \"%(prefix)s.\"$1\".bcf\"}'` -a -d all | bcftools sort -Ob -o %(primer_bcf)s" % self.params
		run_cmd(cmd)
		rm_files(["%s.%s.bcf" % (self.prefix,x) for x in primer_ids])
		rm_files(["%s.%s.bam" % (self.prefix,x) for x in primer_ids])
		rm_files(["%s.%s.bam.bai" % (self.prefix,x) for x in primer_ids])
	def get_calling_params(self):
		dp = []
		cmd = "samtools depth %(bam_file)s" % self.params
		log("Optimising call method")
		for l in subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout:
			arr = l.rstrip().split()
			dp.append(int(arr[2]))
		med_dp = np.median(dp)
		log("Median depth: %s" % med_dp)
		if med_dp<30:
			log("Using low depth approach")
			return "low"
		else:
			log("Using high depth approach")
			return "high"
	def gbcf(self,call_method="optimise",max_dp=None,min_dp=10,threads=4,vtype="snps",bed_file=None,platform="illumina",primers=None,overlap_search=True,chunk_size=50000,mpileup_options=None,low_dp_as_missing=False):
		"""
		Create a gVCF file (for a description see:https://sites.google.com/site/gvcftools/home/about-gvcf)

		Args:
			ref_file(str): reference file (not required if passed to the bam initiator).
			call_method(str): optimise variant calling based on high or low depth. Options: high|low|optimise
			min_dp(int): Minimum depth required to group site into reference-block
		"""
		self.params["min_dp"] = min_dp
		self.params["max_dp"] = max_dp
		self.params["bcf_file"] = "%s.gbcf" % self.prefix
		self.params["bed_file"] = bed_file
		self.params["chunk_size"] = chunk_size
		self.params["cmd_split_chr"] = "splitchr.py %(ref_file)s %(chunk_size)s --bed %(bed_file)s --reformat" % self.params if bed_file else "splitchr.py %(ref_file)s %(chunk_size)s --reformat" % self.params
		self.params["threads"] = threads

		if primers:
			self.params["primer_bed_file"] = "%(prefix)s.primers.bed" % self.params
			TMP = open(self.params["primer_bed_file"],"w")
			positions = self.ref_fa.find_primer_positions(primers)
			for x in sorted(positions,key=lambda d:positions[d]["start"]):
				p = positions[x]
				if p["start"] > p["end"]:
					p["start"],p["end"] = p["end"],p["start"]
				TMP.write("%s\t%s\t%s\t%s\n" % (p["chrom"],p["start"],p["end"],x))
			TMP.close()

		if vtype=="snps": self.params["vtype"] = "-V indels"
		elif vtype=="indels": self.params["vtype"] = "-V snps"
		elif vtype=="both":	self.params["vtype"] = ""
		else: sys.stderr.write("Please provide valid vtype: [snps|indels|both]...Exiting!"); quit(1)
		self.params["primer_cmd"] = " -T ^%(primer_bed_file)s" % self.params if primers else ""

		if call_method=="optimise" and platform=="illumina": call_method = self.get_calling_params()
		self.params["mpileup_options"] = ""
		if platform=="illumina" and  call_method=="high":
			self.params["mpileup_options"] = "-B -a DP,AD"
		elif platform=="illumina" and call_method=="low":
			self.params["mpileup_options"] = "-ABq0 -Q0 -a DP,AD"
		elif platform=="minION":
			if vtype=="snps":
				self.params["mpileup_options"] = "-BIq8 -a DP,AD"
			else:
				self.params["mpileup_options"] = "-Bq8 -a DP,AD"
		else:
			log("Please choose a valid platform...Exiting!",ext=True)
		if mpileup_options:
			self.params["mpileup_options"] = mpileup_options
		self.params["min_dp_cmd"] = "| bcftools filter -e 'FMT/DP<%(min_dp)s' -Ou -S ." % self.params if low_dp_as_missing else ""
		self.params["max_dp_cmd"] = "| bcftools filter -e 'FMT/DP>%(max_dp)s' -Ou -S ." % self.params if max_dp else ""
		cmd = "%(cmd_split_chr)s | parallel --progress --col-sep '\\t' -j %(threads)s \"bcftools mpileup  -f %(ref_file)s %(bam_file)s %(mpileup_options)s -r {1} | bcftools call %(primer_cmd)s %(vtype)s -mg %(min_dp)s | bcftools norm -f %(ref_file)s %(min_dp_cmd)s %(max_dp_cmd)s | bcftools view -Ob -o %(prefix)s_{2}.bcf \"" % self.params
		run_cmd(cmd)
		cmd = "%(cmd_split_chr)s | awk '{print \"%(prefix)s_\"$2\".bcf\"}' | parallel -j  %(threads)s \"bcftools index {}\"" % self.params
		run_cmd(cmd)

		if primers:
			self.params["non_primer_bcf"] = "%(prefix)s.non_primer.bcf" % self.params
			self.params["primer_bcf"] = "%(prefix)s.primer.bcf" % self.params

			if overlap_search:
				#self.params["primer_bam"] = "%(prefix)s.primers.bam" % self.params
				self.generate_primer_bcf()
				#index_bam(self.params["primer_bam"])
				#cmd = "bcftools mpileup  -f %(ref_file)s %(primer_bam)s %(mpileup_options)s -R %(primer_bed_file)s | bcftools call -T %(primer_bed_file)s %(vtype)s -mg %(min_dp)s | bcftools norm -f %(ref_file)s  | bcftools +setGT -Ob -o %(primer_bcf)s -- -t q -i 'FMT/DP<%(min_dp)s' -n ." % self.params
			else:
				cmd = "bcftools mpileup  -f %(ref_file)s %(bam_file)s %(mpileup_options)s -B -R %(primer_bed_file)s | bcftools call %(vtype)s -m | bcftools +setGT -Ob -o %(primer_bcf)s -- -t a -n ." % self.params
				run_cmd(cmd)
			cmd = "bcftools concat -aD -Ob -o %(non_primer_bcf)s `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$2\".bcf\"}'`" % self.params
			run_cmd(cmd)
			cmd = "bcftools concat %(primer_bcf)s %(non_primer_bcf)s | bcftools sort -Ob -o %(bcf_file)s " % self.params
			run_cmd(cmd)
		else:
			cmd = "bcftools concat -aD -Ob -o %(bcf_file)s `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$2\".bcf\"}'`" % self.params
			run_cmd(cmd)

		cmd = "rm `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$2\".bcf*\"}'`" % self.params
		run_cmd(cmd)
		if primers:
			rm_files([self.params["non_primer_bcf"],self.params["primer_bcf"]])

		return bcf(self.params["bcf_file"],prefix=self.prefix)

	def call_variants(self,gff_file=None,bed_file=None,call_method="optimise",min_dp=10,threads=4,mixed_as_missing=False):
		self.params["min_dp"] = min_dp
		self.params["bed_file"] = bed_file
		self.params["cmd_split_chr"] = "splitchr.py %(ref_file)s 50000 --bed %(bed_file)s" % self.params if bed_file else "splitchr.py %(ref_file)s 50000" % self.params
		self.params["gbcf_file"] = "%s.gbcf" % self.prefix
		self.params["missing_bcf_file"] = "%s.missing.bcf" % self.prefix
		self.params["mixed_cmd"] = " bcftools +setGT -- -t q -i 'GT=\"het\"' -n . | bcftools view -e 'F_MISSING==1' |" % self.params if mixed_as_missing else ""
		self.gbcf(call_method=call_method,min_dp=min_dp,threads=threads,vtype="both",bed_file=bed_file,low_dp_as_missing=True)
		self.params["bcf_file"] = "%s.bcf" % self.prefix

		self.params["del_bed"] = bcf(self.params["gbcf_file"]).del_pos2bed()
		view_cmd = "bcftools view %(gbcf_file)s |" % self.params
		mix_cmd = " bcftools +setGT -- -t q -i 'GT=\"het\"' -n . |" % self.params if mixed_as_missing else ""
		out_cmd = "bcftools view -T ^%(del_bed)s -g miss -O b -o %(missing_bcf_file)s" % self.params
		cmd = "%s %s %s" % (view_cmd,mix_cmd,out_cmd)
		run_cmd(cmd)
		out_cmd = "bcftools view -g ^miss -c 1 -O b -o %(bcf_file)s" % self.params
		cmd = "%s %s %s" % (view_cmd,mix_cmd,out_cmd)
		run_cmd(cmd)
		final_bcf = self.params["bcf_file"]
		if gff_file and filecheck(gff_file):
			self.params["gff_file"] = gff_file
			self.params["ann_bcf_file"] = "%(prefix)s.csq.bcf" % self.params
			cmd = "bcftools csq -p m -f %(ref_file)s -g %(gff_file)s %(bcf_file)s -Ob -o %(ann_bcf_file)s" % self.params

			run_cmd(cmd)
			final_bcf = self.params["ann_bcf_file"]


		return bcf(final_bcf,prefix=self.prefix)
	def create_dummy_low_dp_bcf(self,gff_file,min_dp=10,bed_file=None):
		self.params["gff_file"] = gff_file

		contig_line = "\n".join(["##contig=<ID=%s,length=%s>" % (s,len(self.ref_fa_dict[s])) for s in self.ref_fa_dict])
		header = """##fileformat=VCFv4.1
##source=htsbox-pileup-r340
##reference=%s
##contig=<ID=Chromosome,length=4411532>
%s
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
""" % (self.ref_file,contig_line)
		self.params["dp_vcf_file"] = "%(prefix)s.low_cov.vcf" % self.params
		OUT = open(self.params["dp_vcf_file"],"w")
		OUT.write(header)
		self.do_pileup(bed_file)
		for l in open(self.params["temp_pileup"]):
			row = l.rstrip().split()
			ref = row[2]
			alleles = row[3].split(",")
			depth = [int(x) for x in row[4].split(":")[1].split(",")]
			tot_dp = sum(depth)
			if tot_dp>min_dp: continue
			tmp = ["A","C","G","T"]
			fake_allele = tmp.pop()
			if fake_allele==ref: fake_allele = tmp.pop()
			OUT.write("Chromosome\t%s\t.\t%s\t%s\t255\t.\t.\tGT\t1\n" % (row[1],ref,fake_allele))
		OUT.close()
		self.params["ann_bcf_file"] = "%(prefix)s.low_cov.bcf" % self.params
		cmd = "bcftools csq -p m %(dp_vcf_file)s -f %(ref_file)s -g %(gff_file)s -Ob -o %(ann_bcf_file)s" % self.params
		run_cmd(cmd)
		return bcf(self.params["ann_bcf_file"])


	def get_bam_qc(self,cov_thresholds=[1,5,10,20]):
		"""
		Get a qc_bam object

		Args:
			cov_thresholds(list): List of integers to use in the percentage genome covered calculation
		Returns:
			qc_bam: A qc_bam object
		"""
		return qc_bam(self.params["bam_file"],self.params["ref_file"],cov_thresholds)

	def do_pileup(self,bed_file=None):
		self.params["temp"] = bed_file
		self.params["temp_pileup"] = "%(prefix)s.temp.pileup" % self.params
		self.params["temp_bam"] = "%(prefix)s.temp.bam" % self.params
		if bed_file:
			cmd = "htsbox pileup -b %(temp)s -f %(ref_file)s -Q 8 %(bam_file)s > %(temp_pileup)s" % self.params
		else:
			cmd = "htsbox pileup -f %(ref_file)s -Q 8 %(bam_file)s > %(temp_pileup)s" % self.params
		run_cmd(cmd)
	def htsbox_calls(self,bed_file=None):
		self.do_pileup(bed_file=bed_file)
		if bed_file:
			bed_pos = set()
			for l in open(bed_file):
				arr = l.rstrip().split()
				for i in range(int(arr[1]),int(arr[2])+1):
					bed_pos.add((arr[0],str(i)))

		final_calls = defaultdict(lambda :defaultdict(list))
		if self.params["platform"] == "Illumina":
			for l in open(self.params["temp_pileup"]):
				arr = l.rstrip().split()
				if bed_file and (arr[0],arr[1]) not in bed_pos: continue
				calls = arr[3].split(",")
				cov = [int(x) for x in arr[4].split(":")[1].split(",")]
				tot = sum(cov)
				for i in range(len(calls)):
					final_calls[arr[0]][arr[1]].append((calls[i],cov[i]/tot,cov[i]))
		elif self.params["platform"] == "minION":
			for l in open(self.params["temp_pileup"]):
				#Chromosome	  23	  G	   G,G-1C,G+3AAA   0/1:49,1,1
				arr = l.rstrip().split()
				if bed_file and (arr[0],arr[1]) not in bed_pos: continue
				alleles = arr[3].split(",")
				depth = [int(x) for x in arr[4].split(":")[1].split(",")]
				max_allele_dp = max(depth)
				max_allele = alleles[depth.index(max_allele_dp)]
				max_allele_frac = max_allele_dp/tot_dp
				if len(max_allele)>1:
					max_allele = recode_indels([arr[1],max_allele])[1][0]
				if tot_dp<min_dp:
					call = "N"
				if max_allele_frac<min_frac:
					call = "N"
				else:
					call = max_allele
				final_calls[arr[0]][arr[1]].append((call,1,max_allele_dp))
		return final_calls
	def pileup2vcf(self,min_het_frac=0.3,min_hom_frac=0.6,min_dp=10,bed_file=None,indels=True):
		self.params["contig_line"] = "\n".join(["##contig=<ID=%s,length=%s>" % (s,len(self.ref_fa_dict[s])) for s in self.ref_fa_dict])
		header = """##fileformat=VCFv4.1
##source=htsbox-pileup-r340
##reference=%(ref_file)s
%(contig_line)s
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%(prefix)s
""" % self.params
		self.params["temp_pileup"] = "%s.temp.pileup" % self.prefix
		if bed_file:
			self.do_pileup(bed_file)
			bed_pos = set()
			for l in open(bed_file):
				arr = l.rstrip().split()
				for i in range(int(arr[1]),int(arr[2])):
					bed_pos.add(str(i))
		else:
			self.do_pileup()
			pass
		variants = []
		self.params["vcf_file"] = "%s.vcf" % self.prefix
		OUT = open("%(vcf_file)s" % self.params,"w")
		OUT.write(header)
		ref_run_start_pos = -1
		ref_run_start_ref = "X"
		ref_run_min_dp = 0
		for l in open(self.params["temp_pileup"]):
			#Chromosome	  23	  G	   G,G-1C,G+3AAA   0/1:49,1,1
			arr = l.rstrip().split()
			if bed_file:
				if arr[1] not in bed_pos:
					continue
			alleles = arr[3].split(",")
			depth = [int(x) for x in arr[4].split(":")[1].split(",")]
			tot_dp = sum(depth)
			ref = arr[2]
			if ref_run_start_pos==-1:
				ref_run_start_pos = arr[1]
				ref_run_start_ref = ref
				ref_run_min_dp = tot_dp
			max_allele_dp = max(depth)
			max_allele = alleles[depth.index(max_allele_dp)]
			max_allele_frac = max_allele_dp/tot_dp
			adjusted_allele_frac = max_allele_dp/(max_allele_dp+sorted(depth)[-2]) if len(depth)>1 else max_allele_frac
			ref_depth = depth[alleles.index(ref)] if ref in alleles else 0
			if len(max_allele)>1:
				indel = recode_indels([max_allele])
				max_allele = indel[1][0]
				ref = indel[0]

			DP4 = "0,%s,0,%s" % (ref_depth,tot_dp-ref_depth)

			call = max_allele
#			if arr[1]=="13228": import pdb; pdb.set_trace()
			if not indels and (len(max_allele)>1 or len(ref)>1): #INDELS!!!!
				ref_run_end_pos = arr[1]
				if tot_dp<ref_run_min_dp: ref_run_min_dp = tot_dp
			elif tot_dp<min_dp:
				OUT.write("%s\t%s\t.\t%s\t%s\t255\t.\tDP4=%s\tGT:DP\t%s:%s\n" % (arr[0],arr[1],ref,call,DP4,"./.",tot_dp))
				ref_run_start_pos = -1
			elif tot_dp>=min_dp and adjusted_allele_frac>min_hom_frac and call==ref: #REF base call
				ref_run_end_pos = arr[1]
				if tot_dp<ref_run_min_dp: ref_run_min_dp = tot_dp
			elif tot_dp>=min_dp and adjusted_allele_frac<=min_hom_frac and adjusted_allele_frac>min_het_frac: # mixed call
				if call==ref:
					call=alleles[depth.index(sorted(depth)[-2])]
				if (len(call)>1 or len(ref)>1) and not indels:
					ref_run_end_pos = arr[1]
					if tot_dp<ref_run_min_dp: ref_run_min_dp = tot_dp
				else:
					gt="0/1"
					if tot_dp<ref_run_min_dp: ref_run_min_dp = tot_dp
					OUT.write("%s\t%s\t.\t%s\t.\t.\t.\tEND=%s;MinDP=%s\tGT:DP\t0/0:%s\n" % (arr[0],ref_run_start_pos,ref_run_start_ref,int(arr[1])-1,ref_run_min_dp,ref_run_min_dp))
					OUT.write("%s\t%s\t.\t%s\t%s\t255\t.\tDP4=%s\tGT:DP\t%s:%s\n" % (arr[0],arr[1],ref,call,DP4,gt,tot_dp))
					ref_run_start_pos = -1
					variants.append((arr[0],arr[1],ref,call,tot_dp,gt))
			else:
				if call==ref:
					call="."
					gt="0/0"
				else:
					gt="1/1"
					OUT.write("%s\t%s\t.\t%s\t.\t.\t.\tEND=%s;MinDP=%s\tGT:DP\t0/0:%s\n" % (arr[0],ref_run_start_pos,ref_run_start_ref,int(arr[1])-1,ref_run_min_dp,ref_run_min_dp))
				OUT.write("%s\t%s\t.\t%s\t%s\t255\t.\tDP4=%s\tGT:DP\t%s:%s\n" % (arr[0],arr[1],ref,call,DP4,gt,tot_dp))
				ref_run_start_pos = -1
				variants.append((arr[0],arr[1],ref,call,tot_dp,gt))
		OUT.close()
		return variants
	def sambamba_depth(self,outfile,zero_start=False):

		index_bam(self.params["bam_file"])
		fdict = fasta(self.params["ref_file"]).fa_dict
		cov = {}
		self.params["tmp"] = "%s.tmp" % self.prefix
		for s in fdict:
			cov[s] = ["%s\t%s\t0\t0\t0\t0\t0\t0\t0\t%s" % (s,i+1,self.prefix) for i in range(len(fdict[s]))]

		cmd = "sambamba depth base -q 20 -z -t %(threads)s %(bam_file)s > %(tmp)s" % self.params
		run_cmd(cmd)
		for l in open(self.params["tmp"]):
			row = l.rstrip().split()
			if row[0]=="REF": continue
			if zero_start:
				cov[row[0]][int(row[1])] = "\t".join(row)
			else:
				row[1] = str(int(row[1])+1)
				cov[row[0]][int(row[1])-1] = "\t".join(row)

		O = open(outfile,"w")
		for s in fdict:
			for i in range(len(fdict[s])):
					O.write("%s\n" % cov[s][i])
		O.close()
	def get_bed_gt(self,bed_file):
		add_arguments_to_self(self,locals())
		cmd = "bcftools mpileup -f %(ref_file)s -R %(bed_file)s %(bam_file)s -a AD | bcftools call -m | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
		results = defaultdict(lambda : defaultdict(dict))
		for l in cmd_out(cmd):
			#Chromosome	4348079	0/0	51
			chrom,pos,ref,alt,gt,ad = l.rstrip().split()
			pos =int(pos)
			d = {}
			alts = alt.split(",")
			ad = [int(x) for x in ad.split(",")]
			if gt=="0/0":
				d[ref] = ad[0]
			else:
				for i,a in enumerate([ref]+alts):
					d[a] = ad[i]
			results[chrom][pos] = d
		return results
