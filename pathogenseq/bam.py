from __future__ import division
from files import *
from collections import defaultdict
import re
import numpy as np
from qc import *
from utils import *
from mvcf import *
import vcf


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
		index_bam(bam_file,threads=threads)
		if filecheck(bam_file): self.params["bam_file"] = bam_file
		self.params["prefix"] = prefix
		self.prefix = prefix
		if filecheck(ref_file): self.params["ref_file"] = ref_file
		self.params["platform"] = platform
		self.params["threads"] = threads
	def get_calling_params(self):
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
			return "low"
		else:
			print "Using high depth approach"
			return "high"
	def gvcf(self,call_method="optimise",min_dp=10,threads=4):
		"""
		Create a gVCF file (for a description see:https://sites.google.com/site/gvcftools/home/about-gvcf)

		Args:
			ref_file(str): reference file (not required if passed to the bam initiator).
			call_method(str): optimise variant calling based on high or low depth. Options: high|low|optimise
			min_dp(int): Minimum depth required to group site into reference-block
		"""
		self.params["min_dp"] = min_dp
		self.params["vcf_file"] = "%s.gvcf.gz" % self.prefix

		if call_method=="optimise":
			call_method = self.get_calling_params()

		if call_method=="high":
			cmd = "splitchr.py %(ref_file)s 50000 | xargs -i -P %(threads)s sh -c \"samtools mpileup  -ugf %(ref_file)s %(bam_file)s -t DP,AD -r {} | bcftools call -mg %(min_dp)s -V indels -Oz -o %(prefix)s_{}.vcf.gz\"" % self.params
			run_cmd(cmd)
			cmd = "bcftools concat -Oz -o %(vcf_file)s `splitchr.py %(ref_file)s 50000  | awk '{print \"%(prefix)s_\"$1\".vcf.gz\"}'`" % self.params
			run_cmd(cmd)
			cmd = "rm `splitchr.py %(ref_file)s 50000  | awk '{print \"%(prefix)s_\"$1\".vcf.gz\"}'`" % self.params
			run_cmd(cmd)
#			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -t DP | bcftools call -mg %(min_dp)s -V indels -Oz -o %(vcf_file)s" % self.params
		else:
			cmd = "splitchr.py %(ref_file)s 50000 | xargs -i -P %(threads)s sh -c \"samtools mpileup  -ugf %(ref_file)s %(bam_file)s  -aa -ABq0 -Q0 -t DP,AD -r {} | bcftools call -mg %(min_dp)s -V indels -Oz -o %(prefix)s_{}.vcf.gz\"" % self.params
			run_cmd(cmd)
			cmd = "bcftools concat -Oz -o %(vcf_file)s `splitchr.py %(ref_file)s 50000  | awk '{print \"%(prefix)s_\"$1\".vcf.gz\"}'`" % self.params
			run_cmd(cmd)
			cmd = "rm `splitchr.py %(ref_file)s 50000  | awk '{print \"%(prefix)s_\"$1\".vcf.gz\"}'`" % self.params
			run_cmd(cmd)
#			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -ABq0 -Q0 -t DP | bcftools call -mg %(min_dp)s -V indels -Oz -o %(vcf_file)s" % self.params

		variants = []
		vcf_reader = vcf.Reader(open(self.params["vcf_file"]))
		for r in vcf_reader:
			s = r.samples[0]
			variants.append((r.CHROM,r.POS,r.REF,s.gt_bases,s.data.DP,s.data.GT))
		return variants

	def call_variants(self,gff_file=None,bed_file=None,call_method="optimise",min_dp=10,threads=4,mixed_as_missing=False):
		self.params["min_dp"] = min_dp
		self.params["bcf_file"] = "%s.bcf" % self.prefix
		self.params["bed_file"] = bed_file
		self.params["cmd_split_chr"] = "splitchr.py %(ref_file)s 50000 --bed %(bed_file)s" % self.params if bed_file else "splitchr.py %(ref_file)s 50000" % self.params
		self.params["gbcf_file"] = "%s.gbcf" % self.prefix
		self.params["low_dp_bcf_file"] = "%s.low_dp.bcf" % self.prefix

		if call_method=="optimise":
			call_method = self.get_calling_params()

		if call_method=="high":
			cmd = "%(cmd_split_chr)s | parallel -j %(threads)s \"samtools mpileup  -ugf %(ref_file)s %(bam_file)s -B -t DP,AD -r {} | bcftools call -mg %(min_dp)s | bcftools norm -f %(ref_file)s  | bcftools +setGT -Ob -o %(prefix)s_{}.bcf -- -t q -i 'FMT/DP<%(min_dp)s' -n .\"" % self.params
#			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -t DP | bcftools call -mg %(min_dp)s -V indels -Oz -o %(vcf_file)s" % self.params
		else:
			cmd = "%(cmd_split_chr)s | parallel -j %(threads)s \"samtools mpileup  -ugf %(ref_file)s %(bam_file)s  -aa -ABq0 -Q0 -t DP,AD -r {} | bcftools call -mg %(min_dp)s | bcftools norm -f %(ref_file)s | bcftools +setGT -Ob -o %(prefix)s_{}.bcf -- -t q -i 'FMT/DP<%(min_dp)s' -n .\"" % self.params
#			cmd = "samtools mpileup -ugf %(ref_file)s %(bam_file)s -aa -ABq0 -Q0 -t DP | bcftools call -mg %(min_dp)s -V indels -Oz -o %(vcf_file)s" % self.params
		run_cmd(cmd)
		cmd = "%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$1\".bcf\"}' | parallel -j  %(threads)s \"bcftools index {}\"" % self.params
		run_cmd(cmd)
		cmd = "bcftools concat -aD -Ob -o %(gbcf_file)s `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$1\".bcf\"}'`" % self.params
		run_cmd(cmd)
		cmd = "rm `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$1\".bcf*\"}'`" % self.params
		run_cmd(cmd)
		self.params["del_bed"] = bcf(self.params["gbcf_file"]).del_pos2bed()
		cmd = "bcftools view %(gbcf_file)s -T ^%(del_bed)s -g miss -O b -o %(low_dp_bcf_file)s" % self.params
		run_cmd(cmd)
		cmd = "bcftools view %(gbcf_file)s -g ^miss -c 1 -O b -o %(bcf_file)s" % self.params
		run_cmd(cmd)
		final_bcf = self.params["bcf_file"]
		if gff_file and filecheck(gff_file):
			self.params["gff_file"] = gff_file
			self.params["ann_bcf_file"] = "%(prefix)s.csq.vcf.gz" % self.params
			view_cmd = "bcftools view %(bcf_file)s" % self.params
			mixed_cmd = " | bcftools +setGT -- -t q -i 'GT=\"het\"' -n . " % self.params if mixed_as_missing else ""
			csq_cmd = " | bcftools view -e 'F_MISSING==1' |  bcftools csq -p m -f %(ref_file)s -g %(gff_file)s -Ob -o %(ann_bcf_file)s" % self.params

			cmd = "%s %s %s" % (view_cmd,mixed_cmd,csq_cmd)
			run_cmd(cmd)
			final_bcf = self.params["ann_bcf_file"]

		return bcf(final_bcf,prefix=self.prefix)
	def create_dummy_low_dp_bcf(self,gff_file,min_dp=10,bed_file=None):
		self.params["gff_file"] = gff_file
		header = """##fileformat=VCFv4.1
##source=htsbox-pileup-r340
##reference=/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa
##contig=<ID=Chromosome,length=4411532>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""
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
	def pileup2vcf(self,min_het_frac=0.3,min_hom_frac=0.6,min_dp=10,bed_file=None):


		header = """##fileformat=VCFv4.1
##source=htsbox-pileup-r340
##reference=/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa
##contig=<ID=Chromosome,length=4411532>
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
#			if len(max_allele)>1 or len(ref)>1: #INDELS!!!!
#				ref_run_end_pos = arr[1]
#				if tot_dp<ref_run_min_dp: ref_run_min_dp = tot_dp
			if tot_dp>=min_dp and adjusted_allele_frac>min_hom_frac and call==ref:
				ref_run_end_pos = arr[1]
				if tot_dp<ref_run_min_dp: ref_run_min_dp = tot_dp
			elif tot_dp>=min_dp and adjusted_allele_frac<=min_hom_frac and adjusted_allele_frac>min_het_frac:
				if call==ref:
					call=alleles[depth.index(sorted(depth)[-2])]
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
