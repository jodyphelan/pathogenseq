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
# from bokeh.plotting import figure, output_file, show
# from bokeh.layouts import column
from ete3 import Tree
from colour import Color

re_seq = re.compile("([0-9\-]*)([A-Z\*]+)")
re_I = re.compile("([A-Z\*]+)")
number_re = re.compile("[0-9\-]+")

def parse_mutation(x):
	tmp = x.split(">")
	aa_changed = True if len(tmp)>1 else False
	re_obj = re_seq.search(tmp[0])
	change_num = re_obj.group(1)
	ref_aa = re_obj.group(2)
	alt_aa = re_seq.search(tmp[1]).group(2) if aa_changed else None
	return change_num,ref_aa,alt_aa



def load_variants(filename):
	variants = defaultdict(lambda:defaultdict(dict))
	vcf_reader = vcf.Reader(open(filename))
	for rec in tqdm(vcf_reader):
		for s in rec.samples:
			variants[rec.CHROM][rec.POS][s.sample] = s.gt_bases.split("/")[0] if s["GT"]!="./." else "N"
	return variants


def get_missing_positions(bcf_file):
	cmd = "bcftools query -f '%%CHROM\\t%%POS\\n' %s" % bcf_file
	results = []
	for l in subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout:
		row = l.rstrip().split()
		results.append((row[0],row[1]))
	return results






v = True
class bcf:
	def __init__(self,filename,prefix=None,threads=4):
		self.samples = []
		self.filename = filename
		self.threads = threads
		if prefix==None:
			if filename[-4:]==".bcf":
				self.prefix = filename[:-4]
			elif filename[-5:]==".gbcf":
				self.prefix = filename[:-5]
			elif filename[-7:]==".vcf.gz":
				self.prefix = filename[-7:]==".vcf.gz"
			elif filename[-4:]==".vcf":
				self.prefix = filename[-4:]==".vcf"
			else:
				self.prefix = filename
		else:
			self.prefix = prefix
		self.prefix = self.prefix
		self.temp_file = "%s.temp" % self.prefix
		index_bcf(filename,self.threads)
		cmd = "bcftools query -l %(filename)s > %(temp_file)s" % vars(self)
		run_cmd(cmd)
		for l in open(self.temp_file):
			self.samples.append(l.rstrip())
		os.remove(self.temp_file)
		self.vcf = "%s.vcf" % self.prefix

	def del_pos2bed(self):
		self.del_bed = "%s.del_pos.bed" % self.prefix
		OUT = open(self.del_bed,"w")
		cmd = "bcftools view -Ou -v indels %(bcf)s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' | awk 'length($3)>1'" % vars(self)
		sys.stderr.write(cmd)
		j = 0
		for l in subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout:
			j+=1
			row = l.rstrip().split()
			start_pos = int(row[1])+1
			for i in range(start_pos,start_pos+len(row[2])-1):
				OUT.write("%s\t%s\t%s\n" % (row[0],i-1,i))
		if j==0:
			OUT.write("dummy\t1\t1\n")
		OUT.close()

		return self.del_bed



	def load_variants_alt(self,chrom=None,pos=None):
		variants = defaultdict(lambda:defaultdict(dict))
		raw_variants = defaultdict(lambda:defaultdict(dict))
		if chrom and pos:
			cmd = "bcftools view %s %s:%s | bcftools query -f '%%CHROM\\t%%POS[\\t%%IUPACGT]\\n'  | sed 's/\.\/\./N/g'" % (self.filename,chrom,pos)
		else:
			cmd = "bcftools query -f '%%CHROM\\t%%POS[\\t%%IUPACGT]\\n' %s  | sed 's/\.\/\./N/g'" % self.filename
		log(cmd)
		for l in tqdm(subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout):
			row = l.rstrip().split()
			for i in range(len(self.samples)):
				raw_variants[row[0]][row[1]][self.samples[i]] = row[i+2]
		for chrom in raw_variants:
			for pos in raw_variants[chrom]:
				variants[chrom][int(pos)] = raw_variants[chrom][pos]
		if chrom and pos and len(variants)==0:
			log("Variant not found",True)
		if chrom and pos:
			return variants[chrom][int(pos)]
		else:
			return variants

	def load_stats(self):
		self.stats_file = "%s.stats.txt" % self.filename
		cmd = "bcftools stats -v -s - %(filename)s > %(stats_file)s" % vars(self)
		run_cmd(cmd)
		results = defaultdict(lambda:defaultdict(dict))
		for l in open(self.stats_file):
			row = l.rstrip().split("\t")
			if l[0]=="#": continue
			if row[0]=="SN":
				results["SN"][row[2][:-1]] = int(row[3])
			elif row[0]=="AF":
				results["AF"]["SNP"][float(row[2])] = int(row[3])
				results["AF"]["INDEL"][float(row[2])] = int(row[6])
			elif row[0]=="QUAL":
				results["QUAL"]["SNP"][int(row[2])] = int(row[3])
				results["QUAL"]["INDEL"][int(row[2])] = int(row[6])
			elif row[0]=="IDD":
				results["IDD"][int(row[2])] = int(row[3])
			elif row[0]=="ST":
				results["ST"][row[2]] = int(row[3])
			elif row[0]=="DP":
				if row[2][0]==">": continue
				results["DP"][int(row[2])] = int(row[3])
			elif row[0]=="PSC":
				results["PSC"][row[2]]["nRefHom"] = int(row[3])
				results["PSC"][row[2]]["nNonRefHom"] = int(row[4])
				results["PSC"][row[2]]["nHets"] = int(row[5])
		return results
	def plot_stats(self,outfile):
		stats =  self.load_stats()
		output_file(outfile)

		sn = figure(title="Summary stats", x_range=stats["SN"].keys(),toolbar_location=None, tools="")
		sn.vbar(x=stats["SN"].keys(),top=stats["SN"].values(),width=0.9)
		# show the results
		show(sn)

	def split_on_metadata(self,meta_file):
		meta = defaultdict(list)
		for l in open(meta_file):
			#sample	data
			row = l.rstrip().split()
			meta[row[1]].append(row[0])
		for m in meta:
			self.tmp_file = "%s.tmp.txt" % self.prefix
			open(self.tmp_file,"w").write("\n".join(meta[m]))
			self.tmp_bcf = "%s.%s.bcf" % (self.prefix,m)
			cmd = "bcftools view -S %(tmp_file)s %(bcf)s -Ob -o %(tmp_bcf)s" % vars(self)
			run_cmd(cmd)

	def annotate(self,ref_file,gff_file):
		self.ref_file = ref_file
		self.gff_file = gff_file
		self.ann_file = "%s.ann.bcf" % self.prefix
		cmd = "bcftools csq -p m -f %(ref_file)s -g %(gff_file)s %(bcf)s -o %(ann_file)s" % vars(self)
		run_cmd(cmd,verbose=v)

	def extract_matrix(self,matrix_file=None,fmt="old"):
		self.matrix_file = matrix_file if matrix_file==True else self.prefix+".mat"
		if fmt=="new":
			O = open(self.matrix_file,"w").write("chr\tpos\tref\tinfo\ttype\t%s\n" % ("\t".join(self.samples)))
			cmd = "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t.\\t.[\\t%%IUPACGT]\\n' %(bcf)s  | sed 's/\.\/\./N/g' >> %(matrix_file)s" % vars(self)
		elif fmt=="old":
			O = open(self.matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
			cmd = "bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n' %(bcf)s  | sed 's/\.\/\./N/g' >> %(matrix_file)s" % vars(self)
		else:
			log("Choose valid format [old,new]...Exiting!",ext=True)
		run_cmd(cmd,verbose=v)

	def vcf_to_fasta_alt(self,outfile,ref_file,threads=4,chunk_size = 50000, bed_file=None):
		self.ref_file = ref_file
		self.chunk_size = chunk_size
		self.cmd_split_chr = "splitchr.py %(ref_file)s %(chunk_size)s --bed %(bed_file)s --reformat" % vars(self) if bed_file else "splitchr.py %(ref_file)s %(chunk_size)s --reformat" % vars(self)
		self.tmp_file = "%s.tmp.txt" % self.prefix
		self.threads = threads
		cmd = "%(cmd_split_chr)s | parallel --col-sep '\\t' -j %(threads)s \"bcftools view %(filename)s -r {1} -Ou | bcftools query -f '%%POS[\\t%%IUPACGT]\\n' |  datamash transpose > %(prefix)s.{2}.tmp.txt\"" % vars(self)
		run_cmd(cmd)
		cmd = "paste `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$2\".tmp.txt\"}'` > %(tmp_file)s" % vars(self)
		run_cmd(cmd)
		cmd = "rm `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$2\".tmp.txt\"}'`" % vars(self)
		run_cmd(cmd)
		O = open(outfile,"w")
		for i,l in enumerate(open(self.tmp_file)):
			row = l.rstrip().split()
			if i==0: continue
			s = self.samples[i-1]
			seq = "".join(row).replace("./.","N")
			O.write(">%s\n%s\n" % ( s,seq))
		O.close()


	def vcf_to_fasta(self,filename,threads=4):
		"""Create a fasta file from the SNPs"""
		self.threads = threads
		self.tmp_file = "%s.tmp.txt" % self.prefix
		cmd = "bcftools query -f '%%POS[\\t%%IUPACGT]\\n' %(bcf)s |  datamash transpose > %(tmp_file)s" % vars(self)
		run_cmd(cmd)
		O = open(filename,"w")
		for i,l in enumerate(open(self.tmp_file)):
			row = l.rstrip().split()
			if i==0: continue
			s = self.samples[i-1]
			seq = "".join(row).replace("./.","N")
			O.write(">%s\n%s\n" % ( s,seq))
		O.close()

	def bcf2vcf(self):
		if nofile(self.vcf):
			cmd = "bcftools view %(filename)s -Ov -o %(vcf)s" % vars(self)
			run_cmd(cmd)

	def get_variants(self):
		if nofile(self.vcf): self.bcf2vcf()
		vcf_reader = vcf.Reader(open(self.vcf,"r"))
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
		if nofile(self.vcf): self.bcf2vcf()
		vcf_reader = vcf.Reader(open(self.vcf,"r"))
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
		temp_r_script = "%s.temp.R" % self.prefix
		open(temp_r_script,"w").write(rscript)
		cmd = "Rscript %s" % temp_r_script
		run_cmd(cmd)
		rm_files([temp_r_script])

	def merge_in_snps(self,bcf,outfile):
		self.new_bcf = bcf
		self.targets_file = "%(prefix)s.targets" % vars(self)
		self.tmp_file = "%(prefix)s.temp.bcf" % vars(self)
		self.tmp2_file = "%(prefix)s.temp2.bcf" % vars(self)
		self.outfile = outfile
		cmd = "bcftools view -Ou -v snps %(bcf)s | bcftools query -f '%%CHROM\\t%%POS\\n' | awk '{print $1\"\t\"$2-1\"\t\"$2}' > %(targets_file)s" % vars(self)
		run_cmd(cmd)
		cmd = "bcftools view -T %(targets_file)s %(new_bcf)s -Ob -o %(tmp_file)s" % vars(self)
		run_cmd(cmd)
		index_bcf(self.tmp_file,self.threads)
		cmd = "bcftools view -T %(targets_file)s %(bcf)s -Ob -o %(tmp2_file)s" % vars(self)
		run_cmd(cmd)
		index_bcf(self.tmp2_file,self.threads)
		cmd = "bcftools merge --threads %(threads)s -Ou  %(tmp2_file)s %(tmp_file)s | bcftools view -i 'F_MISSING<0.5' -Ob -o %(outfile)s" % vars(self)
		run_cmd(cmd)

	def annotate_from_bed(self,bed_file,outfile=None,nested=False):
		temp_vcf = "%s.temp.vcf" % self.prefix
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
				for i in range(len(tmp)-1):
					if tmp[i] not in tmp[i+1]: switch = False
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
		vcf_reader = vcf.Reader(open(self.vcf))
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
	def bed_subset(self,bed_file,out_file,vcf=False):
		temp_bed = "%s.temp.bed" % self.prefix
		cmd = "awk '{print $1\"\\t\"$2-1\"\\t\"$3}' %s > %s" % (bed_file,temp_bed)
		run_cmd(cmd)
		if vcf:
			cmd = "bcftools view -R %s %s -o %s " % (temp_bed,self.filename,out_file)
		else:
			cmd = "bcftools view -R %s %s -Ob -o %s " % (temp_bed,self.filename,out_file)
		run_cmd(cmd)
		if not vcf:
			return bcf(out_file)
	def odds_ratio(self,bed_file,meta_file,ann_file):
		drugs,meta = load_tsv(meta_file)
		print drugs
		bed_dict = load_bed(bed_file,columns=[5,6],key1=4,key2=5)
		subset_bcf_name = "%s.subset.bcf" % self.prefix
		subset_bcf = self.bed_subset(bed_file,subset_bcf_name)
		variants = subset_bcf.load_csq(ann_file)

		for gene in bed_dict:

			for drug_combo in bed_dict[gene]:
				for var in bed_dict[gene][drug_combo][1].split(";"):
					for drug in bed_dict[gene][drug_combo][0].split(";"):
						if drug not in drugs: continue
						# print("Looking at mutation %s in %s for drug %s" % (var,gene,drug))

						tbl = [[0.5,0.5],[0.5,0.5]]
						change_num,ref_aa,alt_aa = parse_mutation(var)

						if gene not in variants: continue
							# print "%s not in variants"%gene;continue
						if change_num not in variants[gene]:	continue
							# print "%s not in variants[%s]"%(change_num,gene);continue
						try:
							tbl[0][0] += len([s for s in meta.keys() if variants[gene][change_num][s]==alt_aa and meta[s][drug]=="1"])
							tbl[1][0] += len([s for s in meta.keys() if variants[gene][change_num][s]==alt_aa and meta[s][drug]=="0"])
							tbl[0][1] += len([s for s in meta.keys() if variants[gene][change_num][s]==ref_aa and meta[s][drug]=="1"])
							tbl[1][1] += len([s for s in meta.keys() if variants[gene][change_num][s]==ref_aa and meta[s][drug]=="0"])
						except:
							print gene
							print change_num
							print variants[gene][change_num]
							import pdb; pdb.set_trace()
						if tbl[0][0]+tbl[1][0]==1: continue
						OR = (tbl[0][0]/tbl[0][1])/(tbl[1][0]/tbl[1][1])
						print "%s\t%s\t%s\t%s\t%s" % (var,gene,drug,OR,tbl)
#		for record in tqdm(vcf_reader):
#			if record.CHROM in bed_dict and record.POS in bed_dict[record.CHROM]:


	def load_csq(self,ann_file=None,changes=False,use_genomic=True,use_gene=True):
		ann = defaultdict(dict)
		if ann_file:
			for l in tqdm(open(ann_file)):
				#chrom pos gene gene/codon_pos
				row = l.rstrip().split()
				ann[row[0]][int(row[1])] = (row[2],row[3])

		nuc_variants = self.load_variants_alt()
		prot_dict = defaultdict(lambda:defaultdict(dict))
		prot_variants = defaultdict(lambda:defaultdict(dict))
		change_num2pos = defaultdict(lambda:defaultdict(set))
		ref_codons = defaultdict(lambda:defaultdict(dict))

		cmd = "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE\\t%%TBCSQ]\\n' %s" % self.filename
		sys.stderr.write("%s\n"%cmd)
		for line in tqdm(subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout):
			row = line.rstrip().split()
			chrom = row[0]
			pos = int(row[1])
			ref = row[2]
			alt = row[3]
			if chrom in ann and pos in ann[chrom]:
				ann_pos = int(ann[chrom][pos][1])
				ann_gene = ann[chrom][pos][0]
			else:
				ann_pos = None
			if len(row)==4:
				if chrom in ann and pos in ann[chrom]:
					for sample in self.samples:
						prot_variants[ann_gene][ann_pos][sample] = "%s%s>%s" % (ann_pos,ref,alt)
						prot_dict[ann_gene][ann_pos][sample] = nuc_variants[chrom][pos][sample]
						ref_codons[ann_gene][ann_pos] = ref
				continue

			for i in range(4,len(row)-2,3):
				sample = row[i]
				info = row[i+1].split("|") if row[i+1]!="." else row[i+2].split("|")
				if row[i+1][0]=="@": continue
				if info[-1]=="pseudogene": continue
				gene = info[1]
				if info[0]=="intron":continue
				if info[0]=="frameshift&start_lost" or info[0]=="missense&inframe_altering" or info[0]=="missense" or info[0]=="*missense" or info[0]=="start_lost" or info[0]=="*start_lost" or info[0]=="*stop_lost" or info[0]=="stop_lost" or info[0]=="stop_gained" or info[0]=="*stop_gained":
					change_num,ref_aa,alt_aa = parse_mutation(info[5])
					change_num2pos[gene][change_num].add((chrom,pos))
					ref_codons[gene][change_num] = ref_aa
					prot_variants[gene][change_num][row[i]] = info[5]
					prot_dict[gene][change_num][sample] = alt_aa

				elif info[0]=="stop_lost&frameshift" or info[0]=="inframe_insertion" or info[0]=="*inframe_insertion" or info[0]=="inframe_deletion" or info[0]=="*inframe_deletion" or info[0]=="frameshift" or info[0]=="*frameshift" or info[0]=="synonymous" or info[0]=="*synonymous" or info[0]=="stop_retained":
					change_num,ref_nuc,alt_nuc =  parse_mutation(info[6])
					change_num2pos[gene][change_num].add((chrom,pos))
					ref_codons[gene][change_num] = ref_nuc
					change = "%s%s>%s" % (ann_pos,ref_nuc,alt_nuc) if ann_pos else None
					if use_genomic and use_gene and change:
						prot_variants[gene][change_num][row[i]] = change
					elif use_genomic:
						prot_variants[gene][change_num][row[i]] = info[6]
					elif use_gene and change:
						prot_variants[gene][change_num][row[i]] = change
					else:
						prot_variants[gene][change_num][row[i]] = info[5]
					prot_dict[gene][change_num][sample] = alt_nuc
				elif info[0]=="non_coding":
					if chrom in ann and pos in ann[chrom]:
						gene = ann[chrom][pos][0]
						gene_pos = ann[chrom][pos][1]
						prot_variants[gene][gene_pos][sample] = "%s%s>%s" % (gene_pos,ref,alt)
						prot_dict[gene][gene_pos][sample] = alt
						ref_codons[gene][gene_pos] = ref
				else:
					sys.stderr.write(line)
					sys.stderr.write("Unknown variant type...Exiting!\n")
					quit(1)


		for gene in prot_variants:
			for change_num in prot_variants[gene]:
				for s in set(self.samples)-set(prot_variants[gene][change_num].keys()):
					if "N" in  [nuc_variants[chrom][pos][s] for chrom,pos in change_num2pos[gene][change_num]]:
						prot_variants[gene][change_num][s] = "?"
						prot_dict[gene][change_num][s] = "?"
					else:
						pass
						prot_dict[gene][change_num][s] = ref_codons[gene][change_num]
#			if nuc_variants[row[0]][int(row[1])][s]=="N":
#					prot_dict[gene][change_num] = "?"



		for locus in prot_variants:
			prot_variants[locus] = prot_variants[locus].values()
		return prot_variants if changes else prot_dict


	def ancestral_reconstruct(self):
		cmd = "bcftools query -f '%%CHROM\\t%%POS\n' %(bcf)s" % vars(self)
		variants = {}
		for i,l in enumerate(subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout):
			row = l.rstrip().split()
			variants[i] = (row[0],row[1])
		self.reduced_bcf = "%(prefix)s.reduced.bcf" % vars(self)
		cmd = "bcftools view -c 3 %(bcf)s -Ob -o %(reduced_bcf)s" % vars(self)
		run_cmd(cmd)
		reduced = {}
		cmd = "bcftools query -f '%%CHROM\\t%%POS\n' %(reduced_bcf)s" % vars(self)
		for i,l in enumerate(subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout):
			row = l.rstrip().split()
			reduced[i] = (row[0],row[1])
		new_bcf = bcf(self.reduced_bcf)
		self.fasta_file = "%(prefix)s.reduced.snps.fa" % vars(self)
		new_bcf.vcf_to_fasta(self.fasta_file)

		self.tree_file = "%s.newick.txt" % self.prefix
		self.reconstructed_fasta = "%s.reconstructed.fasta" % self.prefix
		cmd = "fastml -s %(fasta_file)s -x %(tree_file)s -j %(reconstructed_fasta)s -qf -mn" % vars(self)
#		run_cmd(cmd,verbose=2)

		fdict = fasta(self.reconstructed_fasta).fa_dict
		t = Tree(self.tree_file, format=1)

		for i in range(len(fdict.values()[0])):
			num_transitions = 0
			for node in t.traverse("postorder"):
				if len(node.get_ancestors())==0: continue
				anc = node.get_ancestors()[0]
				nuc1 = fdict[anc.name][i]
				nuc2 = fdict[node.name][i]
				if nuc1!="?" and nuc2!="?" and nuc1!="N" and nuc2!="N":
					if nuc1!=nuc2:
						num_transitions+=1
						print "%s>%s" % (nuc1,nuc2)
			if num_transitions>1:
				print "Site: %s" % i
				print "Number of transitions: %s" % num_transitions
				print "Location: %s" % (reduced[i][1])
				for node in t.traverse("postorder"):
					nuc = fdict[node.name][i]
					node.add_features(nuc=nuc)
					#p = probs[node.name][i][nuc] if node.name in probs else 1.0
					#node.add_features(prob=p)
				print t.get_ascii(attributes=["name", "nuc"], show_internal=True)

	def itol_from_bcf(self,mutation_file,amino_acid=False):
		if amino_acid:
			all_csq = self.load_csq()
		for l in open(mutation_file):
			mutation = l.rstrip()
			if amino_acid:
				gene,variant = mutation.split("__")
				change_num,ref_aa,alt_aa = parse_mutation(variant)
				if gene in all_csq and change_num in all_csq[gene]:
					variant_dict = all_csq[gene][change_num]
				else:
					continue
			else:
				chrom,pos = mutation.split("__")
				variant_dict = self.load_variants_alt(chrom,pos)

			num_var = len(set(variant_dict.values()))

			cols = [x.get_hex() for x in list(Color("red").range_to(Color("blue"),num_var))]
			col_dict = {d:cols[i] for i,d in enumerate(set(variant_dict.values()))}
			shape_line = "\t".join(["1" for x in range(num_var)])
			col_line = "\t".join(col_dict.values())
			lab_line = "\t".join(col_dict.keys())

			outfile = "%s.itol.txt" % mutation
			OUT = open(outfile,"w")
			OUT.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	%s
COLOR	#ff0000

LEGEND_TITLE	Amino acid
LEGEND_SHAPES	%s
LEGEND_COLORS	%s
LEGEND_LABELS	%s

DATA
""" % (mutation,shape_line,col_line,lab_line))
			for s in self.samples:
				if amino_acid:
					if variant_dict[s]==alt_aa:
						OUT.write("%s\t%s\n" % (s,col_dict[variant_dict[s]]))
				else:
					OUT.write("%s\t%s\n" % (s,col_dict[variant_dict[s]]))

			OUT.close()

	def compress_variants(self):
		cmd = "bcftools query -f '%%CHROM\\t%%POS[\\t%%GT]\\n' %(bcf)s" % vars(self)
		results = defaultdict(list)
		for l in subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout:
			row = l.rstrip().split()
			results[tuple(row[2:])].append([row[0],row[1]])
		final_results = {}
		for i,key in enumerate(results):
			var_name = "variant_%s" % i
			O = open(var_name+".pheno","w")
			for j,s in enumerate(self.samples):
				tmp = "-9"
				if key[j]=="0/0": tmp = "1"
				elif key[j]=="1/1": tmp = "2"
				O.write("0\t%s\t%s\n" % (s,tmp))
			final_results[var_name] = results[key]
			O.close()
		json.dump(final_results,open("%s.compressed_variants.json" % self.prefix,"w"))
		return final_results

	def reheader(self,index_file):
		idx = {}
		for l in open(index_file):
			row = l.rstrip().split()
			if row[0] in idx: sys.stderr.write("Duplicate values in index file (%s)...Exiting!\n"%row[0]); quit(1)
			idx[row[0]] = row[1]

		new_bcf_file = "%(prefix)s.reheader.bcf" % vars(self)
		tmp_header = "%(prefix)s.tmp.header" % vars(self)
		OUT = open(tmp_header,"w")
		for l in subprocess.Popen("bcftools view -h %(bcf)s" % vars(self),shell=True,stdout=subprocess.PIPE).stdout:
			if l[:2]=="##": OUT.write(l); continue
			row = l.rstrip().split()
			for i in range(9,len(row)):
				if row[i] not in idx: print "%s not found in index file...Exiting!" % row[i]; quit(1)
				row[i] = idx[row[i]]
			OUT.write("%s\n" % "\t".join(row))
		OUT.close()
		cmd = "bcftools reheader -h %s %s > %s" % (tmp_header,self.bcf,new_bcf_file)
		run_cmd(cmd)
		rm_files([tmp_header])

	def filt_variants(self,outfile,bed_include=None,bed_exclude=None,threads=4,fmiss=0.1):
		add_arguments_to_self(self,locals())
		print vars(self)
		self.bed_include = "bcftools view -T %s -Ou |" % bed_include if bed_include!=None else ""
		self.bed_exclude = "bcftools view -T ^%s -Ou |" % bed_exclude if bed_exclude!=None else ""
		"""Extract all variant positions"""
		cmd = "bcftools view %(filename)s -Ou | %(bed_include)s %(bed_exclude)s bcftools view --threads %(threads)s -i 'AC>=0 && F_MISSING<%(fmiss)s' -o %(outfile)s -O b" % vars(self)
		run_cmd(cmd)
		return bcf(self.outfile,threads=self.threads)

	def extract_variants(self,outfile,min_dp=10,bed_include=None,bed_exclude=None,threads=4):
		add_arguments_to_self(self,locals())
		self.bed_include = "bcftools view -T %s -Ou |" % bed_include if bed_include!=None else ""
		self.bed_exclude = "bcftools view -T ^%s -Ou |" % bed_exclude if bed_exclude!=None else ""
		cmd = "bcftools +setGT %(filename)s -Ou -- -t q -i 'FMT/DP<%(min_dp)s' -n . | %(bed_include)s %(bed_exclude)s bcftools view --threads %(threads)s -i 'AC>=0' -o %(outfile)s -O b" % vars(self)
		run_cmd(cmd)
		return bcf(self.outfile,threads=self.threads)

	def filt_non_uniq(self,mappability_file,outfile):
		"""Filter out non unique positions"""
		add_arguments_to_self(self,locals())
		non_uniq = []
		self.non_uniq_bed = "%s.genome.non_uniq.bed" % self.prefix
		O = open(self.non_uniq_bed,"w")
		for l in open(self.mappability_file):
			arr = l.rstrip().split()
			if float(arr[3])<1:
				O.write(l)
		O.close()
		cmd = "bcftools view -T ^%(non_uniq_bed)s %(filename)s -O b -o %(outfile)s" % vars(self)
		run_cmd(cmd)
		return bcf(self.outfile,threads=self.threads)

	def sample_filt(self,outfile,miss_cut=0.15,mix_cut=0.15,keep_samples=None):
		"""Filter out low quality samples"""
		add_arguments_to_self(self,locals())
		self.hq_sample_file = "%s.HQ.samples.txt" % self.prefix
		self.lq_sample_file = "%s.LQ.samples.txt" % self.prefix
		self.qual_file = "%s.sample_quals.txt" % self.prefix
		if keep_samples and filecheck(keep_samples):
			self.keep_samples = [x.rstrip() for x in open(keep_samples).readlines()]
		else:
			self.keep_samples = []
		num_calls = int(subprocess.Popen("bcftools view %(filename)s -H | wc -l" % vars(self),shell=True,stdout=subprocess.PIPE).communicate()[0].rstrip())

		miss = {}
		mix = {}
		self.lq_samples = []
		self.hq_samples = []
		HQ = open(self.hq_sample_file,"w")
		LQ = open(self.lq_sample_file,"w")
		QF = open(self.qual_file,"w")
		QF.write("sample\tmix\tmiss\n")

		self.bcftools_stats_file = "%s.bcftools_stats.txt" % self.prefix

		cmd =  "bcftools stats  %(filename)s -s - | grep ^PSC > %(bcftools_stats_file)s" % vars(self)
		run_cmd(cmd)
		for l in open(self.bcftools_stats_file):
			row = l.rstrip().split()
			s = row[2]
			miss[s] = (num_calls-sum([int(row[i]) for i in [3,4,5]]))/num_calls
			mix[s] = int(row[5])/num_calls
			QF.write("%s\t%s\t%s\n" % (s,mix[s],miss[s]))
			if s in self.keep_samples:
				self.hq_samples.append(s)
				HQ.write("%s\n" % s)
			elif miss[s]>self.miss_cut or mix[s]>self.mix_cut:
				self.lq_samples.append(s)
				LQ.write("%s\n" % s)
			else:
				self.hq_samples.append(s)
				HQ.write("%s\n" % s)
		HQ.close()
		LQ.close()
		QF.close()
		cmd = "bcftools view -S %(hq_sample_file)s -a -c 1 -o %(outfile)s -O b %(filename)s" % vars(self)
		run_cmd(cmd)
		return bcf(self.outfile,threads=self.threads)

	def mask_mixed(self,outfile):
		"""Create a BCF file with mixed called masked as missing"""
		add_arguments_to_self(self,locals())
		cmd = "bcftools +setGT %(filename)s -Ou -- -t q -i 'GT=\"het\"' -n . | bcftools view -Ob -o %(outfile)s" % vars(self)
		run_cmd(cmd)
		return bcf(self.outfile,threads=self.threads)
	def generate_consensus(self,ref):
		add_arguments_to_self(self,locals())
		for s in self.samples:
			self.tmp_sample = s
			cmd = "bcftools view -s %(tmp_sample)s -i 'GT==\"./.\"' %(filename)s | bcftools query -f '%%CHROM\\t%%POS\\n'" % vars(self)
			self.tmp_file = "%(prefix)s.%(tmp_sample)s.missing.bed" % vars(self)
			TMP = open(self.tmp_file,"w")
			for l in cmd_out(cmd):
				row = l.rstrip().split()
				TMP.write("%s\t%s\t%s\n" % (row[0],int(row[1])-1,row[1]))
			TMP.close()
			self.tmp_fa = "%(prefix)s.%(tmp_sample)s.tmp.fasta" % vars(self)
			cmd = "bcftools consensus -f %(ref)s %(filename)s -o %(tmp_fa)s -m %(tmp_file)s -s %(tmp_sample)s" % vars(self)
			run_cmd(cmd)
			fa_dict = fasta(self.tmp_fa).fa_dict
			self.final_fa = "%(prefix)s.%(tmp_sample)s.fasta" % vars(self)
			FA = open(self.final_fa,"w")
			for seq in fa_dict:
				FA.write(">%s_%s\n%s" % (self.tmp_sample,seq,fa_dict[seq]))
			FA.close()
			rm_files([self.tmp_file,self.tmp_fa])
