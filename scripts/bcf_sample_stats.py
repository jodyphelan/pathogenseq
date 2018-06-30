#! /usr/bin/env python
import sys
import pathogenseq as ps
import json

infile = sys.argv[1]
ref = sys.argv[2]
outfile = sys.argv[3]
bcf = ps.bcf(infile)
stats = bcf.load_stats(convert=True,ref=ref)
genome_len = sum([len(x) for x in ps.fasta(ref).fa_dict.values()])

print "sample\tnRefHom\tnNonRefHom\tnHets\tnMissing"
for sample in stats["PSC"]:
	s = stats["PSC"][sample]
	s["id"] = sample
	tot_sum = s["nRefHom"]+s["nNonRefHom"]+s["nHets"]
	s["missing"] = genome_len-tot_sum
	print "%s\t%s\t%s\t%s\t%s" % (sample,s["nRefHom"],s["nNonRefHom"],s["nHets"],s["missing"])
	json.dump(s,open(outfile,"w"))
