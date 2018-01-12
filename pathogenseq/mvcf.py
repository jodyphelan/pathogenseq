import sys
from files import *
from fasta import *

v = True
class bcf:
	params = {}
	samples = []
	def __init__(self,filename,prefix=False):
		self.params["bcf"] = filename
		if prefix==False:
			self.params["prefix"] = filename[:-4] if filename[-4:]==".bcf" else filename
		else:
			self.params["prefix"] = prefix
		self.params["temp_file"] = "%s.temp" % self.params["prefix"]
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

	def extract_matrix(self,matrix_file=False):
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
