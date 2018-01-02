import sys
import pathogenseq as ps

if len(sys.argv)!=2:
	print "examl.py <snps.fasta>"
	quit()

snps_file = sys.argv[1]
phylo = ps.phylo(snps_file)
