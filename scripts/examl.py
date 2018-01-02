import sys
import pathogenseq as ps

if len(sys.argv)!=3:
	print "examl.py <snps.fasta> <prefix>"
	quit()

snps_file = sys.argv[1]
prefix = sys.argv[2]
phylo = ps.phylo(snps_file,prefix)
phylo.examl()
