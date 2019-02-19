from __future__ import division
import sys
from .files import *

_chrom_key="chr"
class varmat:
	def __init__(self,filename):
		self.filename = filename
		self.positions = []
		self.header = []
		self.genos = []
		self.samples = []
		self.ref = []
		self.num_samples = None
		self.num_cols = None
		for l in open(self.filename):
			row = l.rstrip().split()
			if row[0]==_chrom_key:
				self.samples = row[3:]
				self.header = row
				self.num_samples = len(row)-3
				self.num_cols = len(row)
				continue
			self.positions.append((row[0],row[1]))
			self.genos.append(row[3:])
			self.ref.append(row[2])
		self.set_positions = set(self.positions)
def compare_varmat(mat1,mat2,report=None):
	if mat1.samples!=mat2.samples:
		log("Samples not identical",ext=True)
	mat1_uniq_pos = []
	mat2_uniq_pos = []
	intersect_pos = []
	geno_discrepancies = {}
	genotype_call_discrepancies = 0
	genotype_position_discrepancies = 0
	for chrompos in mat1.set_positions.union(mat2.set_positions):
		if chrompos in mat1.set_positions and chrompos not in mat2.set_positions:
			mat1_uniq_pos.append(chrompos)
		elif chrompos not in mat1.set_positions and chrompos in mat2.set_positions:
			mat2_uniq_pos.append(chrompos)
		else:
			intersect_pos.append(chrompos)
			mat1_genos = mat1.genos[mat1.positions.index(chrompos)]
			mat2_genos = mat2.genos[mat2.positions.index(chrompos)]
			if mat1_genos!=mat2_genos:
				genotype_position_discrepancies+=1
				tmp = []
				for i in range(mat1.num_samples):
					if mat1_genos[i]!=mat2_genos[i]:
						tmp.append((mat1.samples[i],mat1_genos[i],mat2_genos[i]))
						genotype_call_discrepancies+=1
				geno_discrepancies[chrompos] = tmp
	open(report+".general.txt","w").write("Intersecrion\t%s\nmat1_uniq\t%s\nmat2_uniq\t%s\ngenotype_call_discrepancies\t%s\ngenotype_position_discrepancies\t%s\n" % (len(intersect_pos),len(mat1_uniq_pos),len(mat2_uniq_pos),genotype_call_discrepancies,genotype_position_discrepancies) )
	open(report+".mat1_uniq.txt" ,"w").write("\n".join(["%s\t%s" % (x,y) for x,y in mat1_uniq_pos]))
	open(report+".mat2_uniq.txt" ,"w").write("\n".join(["%s\t%s" % (x,y) for x,y in mat1_uniq_pos]))
	O = open(report+".genotype_discrepancies.txt" ,"w")
	for chrompos in geno_discrepancies:
		for s in geno_discrepancies[chrompos]:
			O.write("%s\t%s\t%s\t%s\t%s\n" % (chrompos[0],chrompos[1],s[0],s[1],s[2]))
	O.close()

	return geno_discrepancies

