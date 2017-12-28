import sys
import subprocess
from files import *
from fasta import *

class phylo:
	params = {}
	def __init__(self,fa_file,prefix,threads=4):
		if filecheck(fa_file):
			self.params["fa_file"] = fa_file
			self.fasta = fasta(fa_file)
			self.params["threads"] = threads
		self.params["prefix"] = prefix
	def examl(self):
		print "Converting to phylip"
		self.params["phylip_file"] = "%s.phylip" % self.params["prefix"]
		self.fasta.write_philip(self.params["phylip_file"])
		print "Creating starting tree"
		cmd = "raxmlHPC-PTHREADS-SSE3 -y -m GTRCAT -p 12345 -s %(fa_file)s  -n StartingTree -T %(threads)s" % self.params
		run_cmd(cmd)
		print "Creating binary file"
		cmd = "parse-examl -s %(phylip_file)s -n %(prefix)s -m DNA" % self.params
		run_cmd(cmd)
		print "Running EXaML"
		cmd = "examl-OMP -s %(prefix)s -n examl -m PSR -D -t RAxML_parsimonyTree.StartingTree" % self.params
