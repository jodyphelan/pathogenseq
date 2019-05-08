from .files import *
from .fasta import *
class delta:
	def __init__(self,filename):
		self.delta = filename
		self.prefix = filename.replace(".delta","")
		self.seqs = {x:fasta(x) for x in open(filename).readline().strip().split()}
	def show_aligns(self,refseq,queryseq):
		add_arguments_to_self(self,locals())
		segments = []
		tmp = []
		for i,l in enumerate(cmd_out("show-aligns %(delta)s %(refseq)s %(queryseq)s -w 105" % vars(self))):
			if l.strip()=="": continue
			row = l.rstrip().split()
			if l[0]==" ":continue
			if i==0:
				seqfiles = row
				continue
			if row[0]=="=========================================================================================================":
				if len(tmp)!=0:
					segments.append(tmp)
					tmp = []
				continue
			tmp.append(row)
		seq_names = (segments[0][3],segments[0][5])
		aln_seqs = {}
		return segments[2:-2:2]
	def get_variants(self):
		add_arguments_to_self(self,locals())
		self.tmp_file = get_random_file()
		run_cmd("show-snps %(delta)s -CTH > %(tmp_file)s" % vars(self))
		variants = []
		for l in open(self.tmp_file):
			row = l.strip().split()
			variants.append({"refseq":row[8],"refpos":int(row[0]),"refnuc":row[1],"queryseq":row[9],"querypos":int(row[3]),"querynuc":row[2]})
		return variants
