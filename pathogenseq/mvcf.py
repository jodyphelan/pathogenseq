from __future__ import division
import sys
import subprocess
from files import *
from fasta import *
import vcf
from collections import defaultdict
import itertools
import json
from tqdm import tqdm

re_ref_aa = re.compile("([0-9]+)([A-Z\*]+)")
re_alt_aa = re.compile("[0-9]+([A-Z\*]+)")


def load_variants(filename):
	variants = defaultdict(lambda:defaultdict(dict))
	vcf_reader = vcf.Reader(open(filename))
	for rec in tqdm(vcf_reader):
		for s in rec.samples:
			variants[rec.CHROM][rec.POS][s.sample] = s.gt_bases.split("/")[0] if s["GT"]!="./." else "N"
	return variants





v = True
class bcf:
	def __init__(self,filename,prefix=None,threads=4):
		self.params = {}
		self.samples = []
		self.params["bcf"] = filename
		self.params["threads"] = threads
		if prefix==None:
			self.params["prefix"] = filename[:-4] if filename[-4:]==".bcf" else filename
		else:
			self.params["prefix"] = prefix
		self.params["temp_file"] = "%s.temp" % self.params["prefix"]
		self.params["vcf"] = "%(prefix)s.vcf" % self.params
		index_bcf(filename,self.params["threads"])
		cmd = "bcftools query -l %(bcf)s > %(temp_file)s" % self.params
		run_cmd(cmd,verbose=v)
		for l in open(self.params["temp_file"]):
			self.samples.append(l.rstrip())
		os.remove(self.params["temp_file"])

	def split_on_metadata(self,meta_file):
		meta = defaultdict(list)
		for l in open(meta_file):
			#sample	data
			row = l.rstrip().split()
			meta[row[1]].append(row[0])
		for m in meta:
			self.params["tmp_file"] = "%s.tmp.txt" % self.params["prefix"]
			open(self.params["tmp_file"],"w").write("\n".join(meta[m]))
			self.params["tmp_bcf"] = "%s.%s.bcf" % (self.params["prefix"],m)
			cmd = "bcftools view -S %(tmp_file)s %(bcf)s -Ob -o %(tmp_bcf)s" % self.params
			run_cmd(cmd)

	def annotate(self,ref_file,gff_file):
		self.params["ref_file"] = ref_file
		self.params["gff_file"] = gff_file
		self.params["ann_file"] = "%s.ann.bcf" % self.params["prefix"]
		cmd = "bcftools csq -f %(ref_file)s -g %(gff_file)s %(bcf)s -o %(ann_file)s" % self.params
		run_cmd(cmd,verbose=v)

	def extract_matrix(self,matrix_file=None):
		self.params["matrix_file"] = matrix_file if matrix_file==True else self.params["prefix"]+".mat"
		O = open(self.params["matrix_file"],"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
		cmd = "bcftools view %(bcf)s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n'  | sed 's/\.\/\./N/g' >> %(matrix_file)s" % self.params
		run_cmd(cmd,verbose=v)

	def vcf_to_fasta_alt(self,filename,threads=4):
		"""Create a fasta file from the SNPs"""
		self.params["threads"] = threads
		cmd = "bcftools query -l %(bcf)s | parallel -j %(threads)s \"(printf '>'{}'\\n' > {}.fa; bcftools view -v snps %(bcf)s | bcftools query -s {} -f '[%%IUPACGT]'  >> {}.fa; printf '\\n' >> {}.fa)\"" % self.params
		run_cmd(cmd)
		O = open(filename,"w")
		for s in self.samples:
			fdict = fasta(s+".fa").fa_dict
			fdict[s] = fdict[s].replace("./.","N")
			O.write(">%s\n%s\n" % ( s,fdict[s]))
			os.remove(s+".fa")
		O.close()

	def vcf_to_fasta(self,filename,threads=4):
		"""Create a fasta file from the SNPs"""
		self.params["threads"] = threads
		self.params["tmp_file"] = "%s.tmp.txt" % self.params["prefix"]
		cmd = "bcftools query -f '%%POS[\\t%%IUPACGT]\\n' %(bcf)s |  datamash transpose > %(tmp_file)s" % self.params
		run_cmd(cmd)
		O = open(filename,"w")
		for i,l in enumerate(open(self.params["tmp_file"])):
			row = l.rstrip().split()
			if i==0: continue
			s = self.samples[i-1]
			seq = "".join(row).replace("./.","N")
			O.write(">%s\n%s\n" % ( s,seq))
		O.close()

	def bcf2vcf(self):
		if nofile(self.params["vcf"]):
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

	def get_venn_diagram_data(self,samples,outfile):
		samples = samples.split(",")
		if len(samples)>4:
			print(samples)
			print("Can't handle more than 4 samples...Exiting!")
			quit()
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
			if si not in self.samples:
				print("Can't find %s in samples...Exiting" % si)
				quit()
			data["id_%s"%i] = si
			data["tot_snps_%s"%i] = tot_snps[si]
		data["outfile"] = outfile
		if len(samples)==2:
			rscript = """
library(VennDiagram)
pdf("%(outfile)s")
draw.pairwise.venn(area1=%(tot_snps_0)s, area2=%(tot_snps_1)s, cross.area=%(overlap_0_1)s, category = c("%(id_0)s","%(id_1)s"),fill=rainbow(2))
dev.off()
""" % data
		elif len(samples)==3:
			rscript = """
library(VennDiagram)
pdf("%(outfile)s")
draw.triple.venn(area1=%(tot_snps_0)s, area2=%(tot_snps_1)s, area3=%(tot_snps_2)s, n12=%(overlap_0_1)s, n23=%(overlap_1_2)s, n13=%(overlap_0_2)s, n123=%(overlap_0_1_2)s, category = c("%(id_0)s","%(id_1)s","%(id_2)s"),fill=rainbow(3))
dev.off()
""" % data
		elif len(samples)==4:
			rscript="""
library(VennDiagram)
pdf("%(outfile)s")
draw.quad.venn(area1=%(tot_snps_0)s, area2=%(tot_snps_1)s, area3=%(tot_snps_2)s, area4=%(tot_snps_3)s,
n12=%(overlap_0_1)s, n13=%(overlap_0_2)s, n14=%(overlap_0_3)s, n23=%(overlap_1_2)s, n24=%(overlap_1_3)s, n34=%(overlap_2_3)s,
n123=%(overlap_0_1_2)s, n124=%(overlap_0_1_3)s, n134=%(overlap_0_2_3)s, n234=%(overlap_1_2_3)s,
n1234=%(overlap_0_1_2_3)s,
category = c("%(id_0)s","%(id_1)s","%(id_2)s","%(id_3)s"),fill=rainbow(4))
dev.off()
""" % data
		temp_r_script = "%s.temp.R" % self.params["prefix"]
		open(temp_r_script,"w").write(rscript)
		cmd = "Rscript %s" % temp_r_script
		run_cmd(cmd)
		rm_files([temp_r_script])
	def merge_in_snps(self,bcf,outfile):
		self.params["new_bcf"] = bcf
		self.params["targets_file"] = "%(prefix)s.targets" % self.params
		self.params["tmp_file"] = "%(prefix)s.temp.bcf" % self.params
		self.params["tmp2_file"] = "%(prefix)s.temp2.bcf" % self.params
		self.params["outfile"] = outfile
		cmd = "bcftools view -v snps %(bcf)s | bcftools query -f '%%CHROM\\t%%POS\\n' | awk '{print $1\"\t\"$2-1\"\t\"$2}' > %(targets_file)s" % self.params
		run_cmd(cmd)
		cmd = "bcftools view -T %(targets_file)s %(new_bcf)s -Ob -o %(tmp_file)s" % self.params
		run_cmd(cmd)
		index_bcf(self.params["tmp_file"],self.params["threads"])
		cmd = "bcftools view -T %(targets_file)s %(bcf)s -Ob -o %(tmp2_file)s" % self.params
		run_cmd(cmd)
		index_bcf(self.params["tmp2_file"],self.params["threads"])
		cmd = "bcftools merge --threads %(threads)s %(tmp2_file)s %(tmp_file)s | bcftools view -i 'F_MISSING<0.5' -Ob -o %(outfile)s" % self.params
		run_cmd(cmd)

	def annotate_from_bed(self,bed_file,outfile=None,nested=False):
		temp_vcf = "%s.temp.vcf" % self.params["prefix"]
		self.vcf_from_bed(bed_file,temp_vcf)
		bed_dict = defaultdict(dict)
		for l in open(bed_file):
			#chrom pos pos allele data
			row = l.rstrip().split()
			bed_dict[row[0]][int(row[1])] = (row[3],row[4])
		vcf_reader = vcf.Reader(open(temp_vcf))
		results = defaultdict(list)
		for record in tqdm(vcf_reader):
			for s in record.samples:
				if s.gt_bases==None: continue
				nuc = s.gt_bases.split("/")[0]
				if nuc==bed_dict[record.CHROM][record.POS][0]:
					results[s.sample].append(bed_dict[record.CHROM][record.POS][1])
		if outfile:
			O = open(outfile,"w")
		for s in self.samples:
			if nested:
				switch = True
				tmp = sorted(list(set(results[s])))
				print tmp
				for i in range(1,len(tmp)):
					if tmp[i] not in tmp[i-1]: switch = False
			else:
				switch = False
			meta = tmp[-1] if switch else ";".join(sorted(list(set(results[s]))))
#			print("%s\t%s" % (s,meta))
			if outfile:
				O.write("%s\t%s\n" % (s,meta))
		if outfile:
			O.close()
		return results
	def extract_compressed_json(self,outfile):
		self.bcf2vcf()
		vcf_reader = vcf.Reader(open(self.params["vcf"]))
		results = defaultdict(lambda: defaultdict(dict))
		for record in tqdm(vcf_reader):
			tmp = defaultdict(list)
			for s in record.samples:
				if s.gt_bases==None:
					tmp["N"].append(self.samples.index(s.sample))
				elif s.gt_nums=="1/1":
					tmp[s.gt_bases.split("/")[0]].append(self.samples.index(s.sample))

			results[record.CHROM][record.POS] = tmp
		json.dump({"variants":results,"samples":self.samples},open(outfile,"w"))
	def vcf_from_bed(self,bed_file,vcf_file):
		temp_bed = "%s.temp.bed" % self.params["prefix"]
		cmd = "awk '{print $1\"\t\"$2-1\"\t\"$3}' %s > %s" % (bed_file,temp_bed)
		run_cmd(cmd)
		cmd = "bcftools view -T %s %s -o %s " % (temp_bed,self.params["bcf"],vcf_file)
		run_cmd(cmd)
	def odds_ratio(self,bed_file,meta_file):
		drugs,meta = load_tsv(meta_file)
		print meta
		self.vcf_from_bed(bed_file,self.params["vcf"])
		bed_dict = load_bed(bed_file,columns=[5,6],key1=4)
		vcf_reader = vcf.Reader(open(self.params["vcf"]))
		variants = self.load_csq()
		for gene in bed_dict:
			for var in bed_dict[gene][1].split(";"):
				for drug in bed_dict[gene][0].split(";"):
					print("Looking at mutation %s in %s for drug %s" % (var,gene,drug))
					if drug not in drugs: continue

					tbl = [[0.5,0.5],[0.5,0.5]]
					codon_num = re_ref_aa.search(var).group(1)
					tbl[0][0] += len([s for s in meta.keys() if s in variants[gene][codon_num] and variants[gene][codon_num][s]==var and meta[s][drug]=="1"])
					tbl[1][0] += len([s for s in meta.keys() if s in variants[gene][codon_num] and variants[gene][codon_num][s]==var and meta[s][drug]=="0"])
					tbl[0][1] += len([s for s in meta.keys() if s not in variants[gene][codon_num] and meta[s][drug]=="1"])
					tbl[1][1] += len([s for s in meta.keys() if s not in variants[gene][codon_num] and meta[s][drug]=="0"])
					print(tbl)
					print((tbl[0][0]/tbl[0][1])/(tbl[1][0]/tbl[1][1]))
#		for record in tqdm(vcf_reader):
#			if record.CHROM in bed_dict and record.POS in bed_dict[record.CHROM]:


	def load_csq(self):

		self.bcf2vcf()
		nuc_variants = load_variants(self.params["vcf"])

		prot_variants = defaultdict(lambda:defaultdict(dict))
		codon_num2pos = defaultdict(lambda:defaultdict(set))
		ref_codons = defaultdict(lambda:defaultdict(dict))
		cmd = "bcftools query -f '%%CHROM\t%%POS[\t%%SAMPLE\t%%TBCSQ]\n' %s" % self.params["vcf"]
		for line in tqdm(subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout):
			row = line.rstrip().split()
			if len(row)==2: continue
			if len(row)==4: continue
			chrom = row[0]
			pos = int(row[1])
			for i in range(2,len(row)-3,3):
				info = row[i+1].split("|")
				if row[i+1][0]=="@": continue
				if info[-1]=="pseudogene": continue
				sample = row[i]
				gene = info[1]
				tmp = info[5].split(">")
				aa_changed = True if len(tmp)>1 else False
				re_obj = re_ref_aa.search(tmp[0])
				codon_num = re_obj.group(1)
				codon_num2pos[gene][codon_num].add((chrom,pos))
				ref_aa = re_obj.group(2)
				ref_codons[gene][codon_num] = ref_aa
				if aa_changed:
					alt_aa = re.search("[0-9]+([A-Z\*]+)",tmp[1]).group(1)
				prot_variants[gene][codon_num][row[i]] = info[5]


		for gene in prot_variants:
			for codon_num in prot_variants[gene]:
				for s in set(self.samples)-set(prot_variants[gene][codon_num].keys()):
					if "N" in  [nuc_variants[chrom][pos][s] for chrom,pos in codon_num2pos[gene][codon_num]]:
						prot_variants[gene][codon_num][s] = "?"
					else:
						#prot_dict[gene][codon_num][s] = ref_codons[gene][codon_num]
						pass


#				print nuc_variants[row[0]][int(row[1])]
#				if nuc_variants[row[0]][row[1]][s]=="N":
#					prot_dict[gene][codon_num] = "?"

#		print prot_dict["Rv0667"]
		return prot_variants
