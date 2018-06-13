#! /usr/bin/env python
import sys
import pathogenseq as ps

args = {}
args["r1"] = sys.argv[1]
args["r2"] = sys.argv[2]
args["centrifuge_db"] = sys.argv[3]
args["fasta_db"] = sys.argv[4]
args["threads"] = sys.argv[5]
args["fq_report"] = ps.get_random_file()
args["log"] = ps.get_random_file()
cmd = "centrifuge -x %(centrifuge_db)s -1 %(r1)s -2 %(r2)s -S %(log)s --report-file %(fq_report)s -p %(threads)s" % args
ps.run_cmd(cmd)

best_ref = ""
best_score = 0
for l in open(args["fq_report"]):
	row = l.rstrip().split("\t")
	if row[4]=="numReads": continue
	if row[2]=="leaf" and int(row[4])>best_score:
		best_ref = row[0]
		best_score = int(row[4])

fasta = ps.fasta(args["fasta_db"]).fa_dict
open("%s.fa" % best_ref,"w").write(">%s\n%s\n" % (best_ref,fasta[best_ref]))
ps.rm_files([args["fq_report"],args["log"]])
