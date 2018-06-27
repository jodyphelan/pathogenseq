from .files import *
from .fasta import *
import multiprocessing

def mdist_worker(filename,split):
	samples = [x.rstrip() for x in cmd_out("bcftools query -l %s" % filename)]
	sample_idx = {s:samples.index(s) for s in samples}
	matrix = [[0 for x in samples] for s in samples]
	miss_matrix = [[0 for x in samples] for s in samples]

	cmd = "bcftools view %s %s |bcftools query -i'GT!=\"ref\"' -f '[\\t%%SAMPLE:%%GT]\\n'" % (filename,split)
	num_snps = 0
	for l in cmd_out(cmd):
		num_snps+=1
		alt_samples = defaultdict(set)
		miss_samples = set()

		row = l.strip().split()
		for x in row:
			s,c = x.split(":")
			if c=="./.":
				miss_samples.add(s)
			else:
				alt_samples[c].add(s)
		for c in alt_samples:
			others  = (set(samples)-alt_samples[c]) - miss_samples
			for s in alt_samples[c]:
				idx = sample_idx[s]
				for x in others:
					matrix[idx][sample_idx[x]]+=1
					matrix[sample_idx[x]][idx]+=1
		for si in miss_samples:
			for sj in set(samples)-miss_samples:
				miss_matrix[sample_idx[si]][sample_idx[sj]]+=1
				miss_matrix[sample_idx[sj]][sample_idx[si]]+=1
			for sj in miss_samples:
				if si==sj: continue
				if sample_idx[si]>sample_idx[sj]:
					miss_matrix[sample_idx[si]][sample_idx[sj]]+=1
	return {"matrix":matrix,"miss_matrix":miss_matrix,"num_snps":num_snps}

def mdist(filename,reffile,outfile,threads=4):
	samples = [x.rstrip() for x in cmd_out("bcftools query -l %s" % filename)]
	sample_idx = {s:samples.index(s) for s in samples}
	matrix = [[0 for x in samples] for s in samples]
	miss_matrix = [[0 for x in samples] for s in samples]
	fa = fasta(reffile)
	jobs = [(filename,x) for x in fa.splitchr(50000,verbose=False)]
	results = []
	with multiprocessing.Pool(threads) as p:
		results = p.starmap(mdist_worker,jobs)
	num_snps = 0
	for r in results:
		num_snps+=r["num_snps"]
		for i in range(len(r["matrix"])):
			for j in range(len(r["matrix"])):
				matrix[i][j] = r["matrix"][i][j]
				matrix[j][i] = r["matrix"][j][i]
				miss_matrix[i][j] = r["miss_matrix"][i][j]
				miss_matrix[j][i] = r["miss_matrix"][j][i]

	for i in range(len(samples)):
		for j in range(len(samples)):
			if j>=i: continue
			scaler = num_snps / (num_snps-miss_matrix[i][j])
			#log("Num SNPs: %s, %s-%s, missing: %s, non missing: %s, abs_dist: %s scale factor: %s, scaled_dist: %s" % (num_snps,self.samples[i],self.samples[j],miss_matrix[i][j],num_snps-miss_matrix[i][j],matrix[i][j],scaler,matrix[i][j]*scaler))

			matrix[i][j] = matrix[i][j]*scaler
			matrix[j][i] = matrix[i][j]
	OUT = open(outfile,"w")
	OUT.write("\t".join(samples)+"\n")
	OUT.write("\n".join(["\t".join([str(d) for d in matrix[j]]) for j in range(len(samples))]))
	OUT.write("\n")
	OUT.close()
	print(matrix)
