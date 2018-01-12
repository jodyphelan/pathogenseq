import sys
import pathogenseq as ps

bam_file = sys.argv[1]
ref_file = sys.argv[2]
chrom = sys.argv[3]
start = sys.argv[4]
end = sys.argv[5]
out_file = sys.argv[6]

start = int(start) if start!="-" else None
end = int(end) if end!="-" else None

x = ps.qc_bam(bam_file,ref_file)
x.plot_cov(chrom,out_file,start=start,end=end)
