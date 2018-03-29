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
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import column
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
		self.params = {}
		self.samples = []
		self.params["bcf"] = filename
		self.bcf = filename
		self.params["threads"] = threads
		if prefix==None:
			self.params["prefix"] = filename[:-4] if filename[-4:]==".bcf" else filename
		else:
			self.params["prefix"] = prefix
		self.prefix = self.params["prefix"]
		self.params["temp_file"] = "%s.temp" % self.params["prefix"]
		self.params["vcf"] = "%(prefix)s.vcf" % self.params
		index_bcf(filename,self.params["threads"])
		cmd = "bcftools query -l %(bcf)s > %(temp_file)s" % self.params
		run_cmd(cmd)
		for l in open(self.params["temp_file"]):
			self.samples.append(l.rstrip())
		os.remove(self.params["temp_file"])

	def del_pos2bed(self):
		self.params["del_bed"] = "%s.del_pos.bed" % self.prefix
		OUT = open(self.params["del_bed"],"w")
		cmd = "bcftools view -v indels %(bcf)s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' | awk 'length($3)>1'" % self.params
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

		return self.params["del_bed"]



	def load_variants_alt(self):
		variants = defaultdict(lambda:defaultdict(dict))
		raw_variants = defaultdict(lambda:defaultdict(dict))
		for l in tqdm(subprocess.Popen("bcftools query -f '%%CHROM\t%%POS[\t%%IUPACGT]\n' %s  | sed 's/\.\/\./N/g'" % self.params["bcf"], shell=True, stdout=subprocess.PIPE).stdout):
			row = l.rstrip().split()
			for i in range(len(self.samples)):
				raw_variants[row[0]][row[1]][self.samples[i]] = row[i+2]
		for chrom in raw_variants:
			for pos in raw_variants[chrom]:
				variants[chrom][int(pos)] = raw_variants[chrom][pos]
		return variants

	def load_stats(self):
		self.params["stats_file"] = "%s.stats.txt" % self.params["bcf"]
		cmd = "bcftools stats -s - %(bcf)s > %(stats_file)s" % self.params
		run_cmd(cmd)
		results = defaultdict(lambda:defaultdict(dict))
		for l in open(self.params["stats_file"]):
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
			self.params["tmp_file"] = "%s.tmp.txt" % self.params["prefix"]
			open(self.params["tmp_file"],"w").write("\n".join(meta[m]))
			self.params["tmp_bcf"] = "%s.%s.bcf" % (self.params["prefix"],m)
			cmd = "bcftools view -S %(tmp_file)s %(bcf)s -Ob -o %(tmp_bcf)s" % self.params
			run_cmd(cmd)

	def annotate(self,ref_file,gff_file):
		self.params["ref_file"] = ref_file
		self.params["gff_file"] = gff_file
		self.params["ann_file"] = "%s.ann.bcf" % self.params["prefix"]
		cmd = "bcftools csq -p m -f %(ref_file)s -g %(gff_file)s %(bcf)s -o %(ann_file)s" % self.params
		run_cmd(cmd,verbose=v)

	def extract_matrix(self,matrix_file=None):
		self.params["matrix_file"] = matrix_file if matrix_file==True else self.params["prefix"]+".mat"
		O = open(self.params["matrix_file"],"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
		cmd = "bcftools view %(bcf)s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n'  | sed 's/\.\/\./N/g' >> %(matrix_file)s" % self.params
		run_cmd(cmd,verbose=v)

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
	def bed_subset(self,bed_file,out_file,vcf=False):
		temp_bed = "%s.temp.bed" % self.params["prefix"]
		cmd = "awk '{print $1\"\\t\"$2-1\"\\t\"$3}' %s > %s" % (bed_file,temp_bed)
		run_cmd(cmd)
		if vcf:
			cmd = "bcftools view -R %s %s -o %s " % (temp_bed,self.params["bcf"],out_file)
		else:
			cmd = "bcftools view -R %s %s -Ob -o %s " % (temp_bed,self.params["bcf"],out_file)
		run_cmd(cmd)
		if not vcf:
			return bcf(out_file)
	def odds_ratio(self,bed_file,meta_file,ann_file):
		drugs,meta = load_tsv(meta_file)
		print drugs
		bed_dict = load_bed(bed_file,columns=[5,6],key1=4,key2=5)
		subset_bcf_name = "%s.subset.bcf" % self.params["prefix"]
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

		cmd = "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE\\t%%TBCSQ]\\n' %s" % self.params["bcf"]
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
		cmd = "bcftools query -f '%%CHROM\\t%%POS\n' %(bcf)s" % self.params
		variants = {}
		for i,l in enumerate(subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout):
			row = l.rstrip().split()
			variants[i] = (row[0],row[1])
		self.params["reduced_bcf"] = "%(prefix)s.reduced.bcf" % self.params
		cmd = "bcftools view -c 3 %(bcf)s -Ob -o %(reduced_bcf)s" % self.params
		run_cmd(cmd)
		reduced = {}
		cmd = "bcftools query -f '%%CHROM\\t%%POS\n' %(reduced_bcf)s" % self.params
		for i,l in enumerate(subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout):
			row = l.rstrip().split()
			reduced[i] = (row[0],row[1])
		new_bcf = bcf(self.params["reduced_bcf"])
		self.params["fasta_file"] = "%(prefix)s.reduced.snps.fa" % self.params
		new_bcf.vcf_to_fasta(self.params["fasta_file"])

		self.params["tree_file"] = "%s.newick.txt" % self.params["prefix"]
		self.params["reconstructed_fasta"] = "%s.reconstructed.fasta" % self.params["prefix"]
		cmd = "fastml -s %(fasta_file)s -x %(tree_file)s -j %(reconstructed_fasta)s -qf -mn" % self.params
#		run_cmd(cmd,verbose=2)

		fdict = fasta(self.params["reconstructed_fasta"]).fa_dict
		t = Tree(self.params["tree_file"], format=1)

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
	def itol_from_bcf(self,mutation_file):
		all_csq = self.load_csq()
		for l in open(mutation_file):
			mutation = l.rstrip()
			gene,variant = mutation.split("__")
			change_num,ref_aa,alt_aa = parse_mutation(variant)
			if gene in all_csq and change_num in all_csq[gene]:
				csq = all_csq[gene][change_num]
			else:
				continue
			num_aa = len(set(csq.values()))

			cols = [x.get_hex() for x in list(Color("red").range_to(Color("blue"),num_aa))]
			col_dict = {d:cols[i] for i,d in enumerate(set(csq.values()))}
			shape_line = "\t".join(["1" for x in range(num_aa)])
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
				if csq[s]==alt_aa:
					OUT.write("%s\t%s\n" % (s,col_dict[csq[s]]))
			OUT.close()

	def compress_variants(self):
		cmd = "bcftools query -f '%%CHROM\\t%%POS[\\t%%GT]\\n' %(bcf)s" % self.params
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
			idx[row[0]] = row[1]
			if row[0] in idx: sys.stderr.write("Duplicate values in index file...Exiting!"); quit(1)
		new_bcf_file = "%(prefix)s.reheader.bcf" % self.params
		tmp_header = "%(prefix)s.tmp.header" % self.params
		OUT = open(tmp_header,"w")
		for l in subprocess.Popen("bcftools view -h %(bcf)s" % self.params,shell=True,stdout=subprocess.PIPE).stdout:
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
