import sys
import pathogenseq as ps
infile = sys.argv[1]
outfile = sys.argv[2]
threads = sys.argv[3]

bcf = ps.bcf(infile)
bcf.vcf_to_fasta(outfile,threads)
