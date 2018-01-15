import sys
import subprocess
from files import *
from fasta import *
import vcf
from collections import defaultdict
import itertools

v = True
class bcf:
	params = {}
	samples = []
	def __init__(self,filename,prefix=None):
		self.params["bcf"] = filename
		if prefix==None:
			self.params["prefix"] = filename[:-4] if filename[-4:]==".bcf" else filename
		else:
			self.params["prefix"] = prefix
		self.params["temp_file"] = "%s.temp" % self.params["prefix"]
		self.params["vcf"] = "%(prefix)s.vcf" % self.params
		cmd = "bcftools query -l %(bcf)s > %(temp_file)s" % self.params
		run_cmd(cmd,verbose=v)
		for l in open(self.params["temp_file"]):
			self.samples.append(l.rstrip())
		os.remove(self.params["temp_file"])

	def annotate(self,ref_file,gff_file):
		self.params["ref_file"] = ref_file
		self.params["gff_file"] = gff_file
		self.params["ann_file"] = "%s.ann.bcf" % self.params["prefix"]
		cmd = "bcftools csq -f %(ref_file)s -g %(gff_file)s %(bcf)s -o %(ann_file)s" % self.params
		run_cmd(cmd,verbose=v)

	def extract_matrix(self,matrix_file=None):
		self.params["matrix_file"] = matrix_file if matrix_file==True else self.params["prefix"]+".mat"
		O = open(self.params["matrix_file"],"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
		cmd = "bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n' %(bcf)s | sed 's/\.\/\./N/g' >> %(matrix_file)s" % self.params
		run_cmd(cmd,verbose=v)
	def vcf_to_fasta(self,filename,threads=20):
		"""Create a fasta file from the SNPs"""
		self.params["threads"] = threads
		cmd = "bcftools query -l %(bcf)s | parallel -j %(threads)s \"(printf '>'{}'\\n' > {}.fa; bcftools query -s {} -f '[%%IUPACGT]' %(bcf)s >> {}.fa; printf '\\n' >> {}.fa)\"" % self.params
		run_cmd(cmd)
		O = open(filename,"w")
		for s in self.samples:
			fdict = fasta(s+".fa").fa_dict
			fdict[s] = fdict[s].replace("./.","N")
			O.write(">%s\n%s\n" % ( s,fdict[s]))
			os.remove(s+".fa")
		O.close()
	def bcf2vcf(self):
		cmd = "bcftools view %(bcf)s -Ov -o %(vcf)s" % self.params
		run_cmd(cmd)
	def get_variants(self):
		if nofile(self.params["vcf"]): self.bcf2vcf()
		vcf_reader = vcf.Reader(open(self.params["vcf"],"r"))
		results = defaultdict(lambda:defaultdict(dict))
		for record in vcf_reader:
			for s in record.samples:
				results[record.CHROM][record.POS][s.sample] = s.gt_bases.split("/")[0]
		return results
	def get_venn_diagram_data(self,samples):
		samples = samples.split(",")
		if nofile(self.params["vcf"]): self.bcf2vcf()
		vcf_reader = vcf.Reader(open(self.params["vcf"],"r"))
		results = defaultdict(int)
		tot_snps = defaultdict(int)
		data = defaultdict(int)
		for record in vcf_reader:
			tmp = []
			for s in record.samples:
				if s.sample not in samples: continue
				if s.gt_nums=="1/1":
					tmp.append(s.sample)
					tot_snps[s.sample]+=1
			for x in itertools.combinations(tmp,2):
				tmp_str = "_".join(sorted([str(samples.index(d)) for d in x]))
				data["overlap_"+tmp_str] +=1
			for x in itertools.combinations(tmp,3):
				tmp_str = "_".join(sorted([str(samples.index(d)) for d in x]))
				data["overlap_"+tmp_str] +=1
			for x in itertools.combinations(tmp,4):
				tmp_str = "_".join(sorted([str(samples.index(d)) for d in x]))
				data["overlap_"+tmp_str] += 1

		for i,si in enumerate(samples):
			data["id_%s"%i] = si
			data["tot_snps_%s"%i] = tot_snps[si]
		if len(samples)==2:
			rscript = """
library(VennDiagram)
draw.pairwise.venn(area1=%(tot_snps_0)s, area2=%(tot_snps_1)s, cross.area=%(overlap_0_1)s, category = c("%(id_0)s","%(id_1)s"),fill=rainbow(2))
""" % data
		elif len(samples)==3:
			rscript = """
library(VennDiagram)
draw.triple.venn(area1=%(tot_snps_0)s, area2=%(tot_snps_1)s, area3=%(tot_snps_2)s, n12=%(overlap_0_1)s, n23=%(overlap_1_2)s, n13=%(overlap_0_2)s, n123=%(overlap_0_1_2)s, category = c("%(id_0)s","%(id_1)s","%(id_2)s"),fill=rainbow(3))
""" % data
		elif len(samples)==4:
			rscript="""
library(VennDiagram)
draw.quad.venn(area1=%(tot_snps_0)s, area2=%(tot_snps_1)s, area3=%(tot_snps_2)s, area4=%(tot_snps_3)s,
n12=%(overlap_0_1)s, n13=%(overlap_0_2)s, n14=%(overlap_0_3)s, n23=%(overlap_1_2)s, n24=%(overlap_1_3)s, n34=%(overlap_2_3)s,
n123=%(overlap_0_1_2)s, n124=%(overlap_0_1_3)s, n134=%(overlap_0_2_3)s, n234=%(overlap_1_2_3)s,
n1234=%(overlap_0_1_2_3)s,
category = c("%(id_0)s","%(id_1)s","%(id_2)s","%(id_3)s"),fill=rainbow(4))

""" % data
		print rscript
#
