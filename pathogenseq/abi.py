from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt
from .files import *
from .nucmer import *
from .mvcf import *
from .fasta import *
import matplotlib as mpl
mpl.use('TkAgg')
class abi:
	def __init__(self,in_obj):
		if isinstance(in_obj,str):
			self.filename = in_obj
			self.record = SeqIO.read(in_obj,'abi')
			self.prefix = in_obj.replace(".ab1","")
		elif isinstance(in_obj,SeqIO.SeqRecord):
			self.record = in_obj
			self.prefix = in_obj.prefix
		self.quals = [ord(x) for x in self.record.annotations["abif_raw"]["PCON1"]]
	def trim_seq(self,rec):
		return abi(SeqIO.AbiIO._abi_trim(self.record))
	def plot_chromatogram(self,start,end,refsequence):
		self.signals = {}
		self.channels = {'DATA9':"A", 'DATA10':"C", 'DATA11':"G", 'DATA12':"T"}
		for c in self.channels:
			self.signals[self.channels[c]] = self.record.annotations['abif_raw'][c]
		self.ploc1 = self.record.annotations["abif_raw"]["PLOC1"]
		cols = {"A":"gold","C":"red","G":"green","T":"blue"}
		signal_start = self.ploc1[start-1]-3
		signal_end = self.ploc1[end]+5
		xvals = list(range(signal_start,signal_end))
		for c in ["A","C","G","T"]:
			plt.plot(xvals,self.signals[c][signal_start:signal_end], color=cols[c])
		plt.xticks(self.ploc1[start-1:end],str(self.record.seq)[start-1:end])
		for i,x in enumerate(refsequence):
			col = "red" if i==4 else "black"
			plt.text(self.ploc1[start-1:end][i],(plt.axis()[3] - plt.axis()[2])*-0.15 ,x,horizontalalignment='center',color=col)
#		plt.fill_between([self.ploc1[start-1],self.ploc1[end]],plt.axis()[2],plt.axis()[3])
		plt.show()
	def write_seq(self):
		self.fasta = "%s.fasta" % self.prefix
		open(self.fasta,"w").write(">%s\n%s\n" % (self.prefix,self.record.seq))

	def nucmer_align(self,refseq):
		add_arguments_to_self(self,locals())
		self.write_seq()
		run_cmd("nucmer %(refseq)s %(fasta)s -p %(prefix)s" % vars(self))
		return delta("%s.delta" % self.prefix)
	def get_variants_vcf(self,refseq):
		add_arguments_to_self(self,locals())
		self.write_seq()
		fa = fasta(self.prefix+".fasta")
		return bcf(fa.get_ref_variants(refseq,self.prefix))
	def get_variants(self,refseq):
		add_arguments_to_self(self,locals())
		variants = []
		for l in cmd_out("minimap2 %(refseq)s %(prefix)s.fasta --cs | sort -k6,6 -k8,8n | paftools.js call -l 100 -L 100 -" % vars(self)):
			row = l.strip().split()
			if row[0]!="V": continue
			variants.append({"refseq":row[1],"refpos":int(row[2]),"refnuc":row[6],"queryseq":row[8],"querypos":int(row[9]),"querynuc":row[7]})
		return variants
	def load_maf(self,refseq):
		self.maf = {}
		for l in cmd_out("minimap2 %(refseq)s %(prefix)s.fasta --cs=long | sort -k6,6 -k8,8n | paftools.js view -f maf -" % vars(self)):
			row = l.strip().split()
			if l=="": continue
			if row[0]!="s": continue
			self.maf[row[1]] = {"start":int(row[2]),"seq":row[6].upper()}
	def get_maf_refseq(self,start,end,refseq = None):
		if "maf" not in vars(self): self.load_maf(refseq)
		start = start - self.maf[self.prefix]["start"]-1
		end = end - self.maf[self.prefix]["start"]
		print((start,end))
		return self.maf["KY241726"]["seq"][start:end]
	def plot_variants(self,refseq):
		variants = self.get_variants(refseq)
		for var in variants:
			print(var)
			refsequence = self.get_maf_refseq(var["querypos"]-3,var["querypos"]+5,refseq)
			self.plot_chromatogram(var["querypos"]-3,var["querypos"]+5,refsequence)
