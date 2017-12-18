import sys
import pathogenseq

bam_file = sys.argv[1]
ref_file = sys.argv[2]
out_file = sys.argv[3]

x = qc_bam(bam_file,ref_file)
x.plot_cov("Chromosome",out_file)
x.save_cov("test.json")
