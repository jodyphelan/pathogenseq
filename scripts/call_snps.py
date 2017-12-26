import pathogenseq as ps

bam_file = sys.argv[1]
ref_file = sys.argv[2]
prefix = sys.argv[3]

bam = ps.bam(bam_file=bam_file,ref_file=ref_file,prefix=prefix)
bam.call_snps()
