import sys
import re
from files import *
from utils import *

class fasta:
	"""
	Class to represent fasta seuqnces in a python dict.

	Args:
		filename(str): Location of the fasta file

	Returns:
		fasta: A fasta class object
	"""
	def __init__(self,filename):
		fa_dict = {}
		seq_name = ""
		self.fa_file = filename
		for l in open(filename):
			line = l.rstrip()
			if line[0] == ">":
				seq_name = line[1:].split()[0]
				fa_dict[seq_name] = []
			else:
				fa_dict[seq_name].append(line)
		result = {}
		for seq in fa_dict:
			result[seq] = "".join(fa_dict[seq])
			result[seq] = result[seq].upper()
		self.fa_dict = result
	def n50(self):
		"""Return n50"""
		numlist = [len(self.fa_dict[x]) for x in self.fa_dict.keys()]
		numlist.sort()
		newlist = []
		for x in numlist :
			newlist += [x]*x
		# take the mean of the two middle elements if there are an even number
		# of elements.  otherwise, take the middle element
		if len(newlist) % 2 == 0:
			medianpos = len(newlist)/2
			return float(newlist[medianpos] + newlist[medianpos-1]) /2
		else:
			medianpos = len(newlist)/2
		return int(newlist[medianpos])
	def seq_names(self):
		"""Return sequence names"""
		return self.fa_dict.keys()
	def total_len(self):
		"""Return total length of all sequences"""
		length = 0
		for s in self.fa_dict:
			length+=len(self.fa_dict[s])
		return length
	def locate_motifs(self,motif,strand="both"):
		"""Get locations of motifs"""
		IUPAC_codes = {"A":"A","C":"C","G":"G","T":"T","R":"[AG]","Y":"[CT]","S":"[GC]","W":"[AT]","K":"[GT]","M":"[AC]","B":"[CGT]","D":"[AGT]","H":"[ACT]","V":"[ACG]","N":"[ACGT]"}
		pattern = "".join([IUPAC_codes[x] for x in motif])
		re_obj = re.compile(pattern)
		result = []
		motif_len = len(motif)
		for s in self.seq_names():
			if strand=="both" or strand=="+":
				res = [i for i in range(len(self.fa_dict[s])-motif_len) if re_obj.match(self.fa_dict[s][i:i+motif_len])!=None]
				for x in res:
					result.append((s,x,"+"))
			if strand=="both" or strand=="-":
				res = [i for i in range(len(self.fa_dict[s])-motif_len) if re_obj.match(revcom(self.fa_dict[s][i:i+motif_len]))!=None]
				for x in res:
					result.append((s,x,"-"))
		return result
	def write_philip(self,out_file):
		fa_dict = self.fa_dict
		O = open(out_file,"w")
		O.write("\t%s\t%s\n" % (len(fa_dict.keys()),len(fa_dict[fa_dict.keys()[0]])))
		for s in fa_dict:
			O.write("%s\t%s\n" % (s,fa_dict[s]))
		O.close()
	def get_seq(self,chrom,start,end=None):
		if end:
			return self.fa_dict[chrom][start-1:end-1]
		else:
			return self.fa_dict[chrom][start-1]
	def get_VCF(self,ref_file,prefix):
		ref_dict = fasta(ref_file)
		params = {"fa_file":self.fa_file,"ref_file":ref_file,"prefix":prefix}
		vcf_file = "%s.vcf" % prefix
		cmd = "dnadiff %(ref_file)s %(fa_file)s -p %(prefix)s" % params
		run_cmd(cmd)
		snps_file = "%s.snps" % prefix
		del_lines = []
		ins_lines = []
		indel_line_set = set()
		tmp = []
		prev_type = None
		prev_pos = None
		prev_chrom = None
		lines = [l.rstrip().split() for l in open(snps_file).readlines()]
		for i in range(len(lines)):
			#4215484	G	C	4214371	2073	196049	4411532	4410422	1	1	Chromosome	WBB445_ARS7496|quiver
			#4212835	C	.	4211722	1	198698	4411532	4410422	1	1	Chromosome	WBB445_ARS7496|quiver
			row = lines[i]
			pos = int(row[0])
			chrom = row[10]
			if row[1]=="." or row[2]==".":
				indel_line_set.add(i)
				if prev_pos==None:
					tmp = [i]
				elif row[1]==".":
					if pos==prev_pos:
						tmp.append(i)
					else:
						if prev_type=="ins":
							ins_lines.append(tmp)
						else:
							del_lines.append(tmp)
						tmp = [i]
				elif  row[2]==".":
					if pos==prev_pos+1:
						tmp.append(i)
					else:
						if prev_type=="ins":
							ins_lines.append(tmp)
						else:
							del_lines.append(tmp)
						tmp = [i]
				prev_type = "ins" if row[1]=="." else "del"
				prev_pos = pos
				prev_chrom = chrom

		variants = defaultdict(dict)

		for indel_pos in del_lines:
			positions = [lines[i][0] for i in indel_pos]
			first_pos = int(positions[0])
			chrom = lines[0][10]
			bases = [lines[i][1] for i in indel_pos]
			switch = True
			start_pos = first_pos-1
			while switch:
				n = ref_dict.get_seq(chrom,start_pos)
				end_base = bases[-1]
				if n==end_base:
					start_pos-=1
					bases.insert(0,bases.pop())
				else:
					switch=False
			alt_seq = ref_dict.get_seq(chrom,start_pos)
			ref_seq = alt_seq+"".join(bases)
			variants[chrom][start_pos] = (ref_seq,alt_seq)
		for indel_pos in ins_lines:
			positions = [lines[i][0] for i in indel_pos]
			first_pos = int(positions[0])
			chrom = lines[0][10]
			bases = [lines[i][2] for i in indel_pos]
			switch = True
			start_pos = first_pos
			while switch:
				n = ref_dict.get_seq(chrom,start_pos)
				end_base = bases[-1]
				if n==end_base:
					start_pos-=1
					bases.insert(0,bases.pop())
				else:
					switch=False
			ref_seq = ref_dict.get_seq(chrom,start_pos)
			alt_seq = ref_seq+"".join(bases)
			variants[chrom][start_pos] = (ref_seq,alt_seq)

		for i in set(range(len(lines)))-indel_line_set:
			row = lines[i]
			pos,ref_seq,alt_seq = row[:3]
			chrom = row[10]
			variants[chrom][int(pos)] = (ref_seq,alt_seq)
		OUT = open(vcf_file,"w")
		OUT.write("""##fileformat=VCFv4.1
##reference=/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa
##contig=<ID=Chromosome,length=4411532>
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s
""" % prefix)
		npos = None
		nchrom = None
		for chrom in sorted(variants):
			for pos in sorted(variants[chrom]):
				if nchrom!=chrom:
					if pos!=1:
						npos=1
						OUT.write("%s\t%s\t.\t%s\t.\t.\t.\tEND=%s;MinDP=20\tGT:DP\t0/0:20\n" % (chrom,npos,ref_dict.get_seq(chrom,npos),pos-1))
					else:
						pass
				else:
					OUT.write("%s\t%s\t.\t%s\t.\t.\t.\tEND=%s;MinDP=20\tGT:DP\t0/0:20\n" % (nchrom,npos+1,ref_dict.get_seq(chrom,npos),pos-1))
				var = variants[chrom][pos]
				OUT.write("%s\t%s\t.\t%s\t%s\t255\t.\t.\tGT:DP\t%s:%s\n" % (chrom,pos,var[0],var[1],1/1,20))
				npos = pos
				nchrom = chrom
		OUT.close()

	def splitchr(self,size):
		lengths = {s:len(self.fa_dict[s]) for s in self.fa_dict}
		for s in self.fa_dict:
			start = 0
			end = start+size
			while end<lengths[s]:
				print "%s:%s-%s" % (s,start+1,end)
				start+=size
				end+=size
			print "%s:%s-%s" % (s,start,lengths[s])


def revcom(s):
	"""Return reverse complement of a sequence"""
	def complement(s):
			basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
			letters = list(s)
			letters = [basecomplement[base] for base in letters]
			return ''.join(letters)
	return complement(s[::-1])
