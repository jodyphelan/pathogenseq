from __future__ import division
import subprocess
from .files import *
import os
from tqdm import tqdm
import random
r = random.SystemRandom()
class ann:
    annfile = ""
    tabix = ""
    def __init__(self,filename,tabix):
        self.annfile = filename
        self.tabix = tabix
    def pos2ann(self,pos_tuple_list):
        #pos_tuple_list = [("Chromosome",1),("Chromosome",2)]
        if len(pos_tuple_list)<5000:
            p1offset = -1
            p2offset = 0
            stype = "R"
        else:
            p1offset = 0
            p2offset = 1
            stype = "T"

        num = r.randint(1,1000000)+r.randint(1,1000000)
        temp_bed_file = "temp.%s.bed" % (num)
        OUT = open(temp_bed_file,"w")
        for chrom,pos in pos_tuple_list:
            OUT.write("%s\t%s\t%s\n" % (chrom,int(pos)+p1offset,int(pos)+p2offset))
        OUT.close()
        results = defaultdict(dict)
        for l in subprocess.Popen("%s %s -%s %s " % (self.tabix,self.annfile,stype,temp_bed_file),stdout=subprocess.PIPE,shell=True).stdout:
            # added .decode() to convert binary output of subprocess.Popen to string
            l = l.decode()
            arr = l.rstrip().split()
            results[arr[0]][int(arr[1])] = {"change_pos":arr[7],"ref_nt":arr[2],"ref_codon":arr[6],"ref_aa":arr[11],"chr":arr[0],"pos":int(arr[1]),"rv":arr[15],"gene":arr[16],"gene_syn":arr[17],"ncr":arr[18],"start":arr[19],"end":arr[20],"strand":arr[21],"codon_num":arr[24],"gene_nt":arr[25],"operon":arr[26]}

        os.remove(temp_bed_file)
        return results


class gemma_results:
	def __init__(self,filename=None):
		self.data = []
		if filename:
			for l in tqdm(open(filename)):
				row = l.rstrip().split()
				if row[1]=="rs":continue
				chrom,pos,alt = row[1].split("_")
				self.data.append({"chrom":chrom,"pos":int(pos),"ref":row[4],"alt":row[5],"pval":float(row[11])})
			self.data = sorted(self.data,key=lambda x:x["pval"])
	def add_results(self,data):
		for d in data:
			self.data.append(d)
		self.data = sorted(self.data,key=lambda x:x["pval"])
	def top_n_hits(self,n=10):
		tmp = gemma_results()
		tmp.add_results(self.data[:n])
		return tmp
	def cutoff_hits(self,cutoff=1e-5):
		tmp = gemma_results()
		tmp.add_results([x for x in self.data if x["pval"]<cutoff])
		return tmp
	def __str__(self):
		print(self.data)
	def tb_annotate_hits(self,ann_file):
		ann_obj = ann(ann_file,"tabix")
		annotation = ann_obj.pos2ann(sorted([(row["chrom"],row["pos"]) for row in self.data],key=lambda x:x[1]))
		for d in self.data:
			d["gene"] = annotation[d["chrom"]][d["pos"]]["rv"]
			d["annotation"] = annotation[d["chrom"]][d["pos"]]
	def restrict_hits_to_gene(self,genes):
		if "gene" not in self.data[0]:
			print("Please annotate hits with tb_annotate_hits()")
			return
		tmp = gemma_results()
		tmp.add_results([x for x in self.data if x["gene"] in genes])
		return tmp
	def create_tb_panel(self,drug,ann_file,cutoff=1e-5):
		target_genes = {
		"amikacin":["rrs"],
		"bedaquiline":["Rv0678"],
		"capreomycin":["rrs","Rv1694"],
		"clofazimine":["Rv0678"],
		"ethambutol":["Rv3794","Rv3795","Rv3793","Rv1267c","Rv3793-Rv3794"],
		"ethionamide":["Rv3854c","Rv3855","Rv1482c-Rv1483","Rv1484"],
		"fluoroquinolones":["Rv0005","Rv0006"],
		"isoniazid":["Rv2428","Rv2427A-Rv2428","Rv1482c-Rv1483","Rv1484","Rv2245","Rv1908c","Rv1908c-Rv1909c"],
		"kanamycin":["Rv2416c-Rv2417c","rss"],
		"linezolid":["Rv0701","rrl"],
		"para-aminosalicylic_acid":["Rv2447c","Rv2671","Rv2764c"],
		"pyrazinamide":["Rv3601c","Rv2043c","Rv2043c-Rv2044c","Rv1630"],
		"rifampicin":["Rv0667","Rv0668"],
		"streptomycin":["Rv3919c","Rv0682","rrs"]}
		self = self.cutoff_hits(cutoff)
		self.tb_annotate_hits(ann_file)
		self = self.restrict_hits_to_gene(target_genes[drug.lower()])
		for d in self.data:
			mut = ""
			if
			print("%s\t%s\t%s\t%s\t%s\t%s" % (drug,d["pos"],d["ref"],d["alt"],d["gene"]))
