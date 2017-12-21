import sys
import re

class fasta:
	"""
	Class to represent fasta seuqnces in a python dict.

	Args:
		filename(str): Location of the fasta file

	Returns:
		fasta: A fasta class object
	"""
	fa_dict = {}
	def __init__(self,filename):
		fa_dict = {}
		seq_name = ""
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

def revcom(s):
	"""Return reverse complement of a sequence"""
	def complement(s):
			basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
			letters = list(s)
			letters = [basecomplement[base] for base in letters]
			return ''.join(letters)
	return complement(s[::-1])
