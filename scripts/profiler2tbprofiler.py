#! /usr/bin/env python
import sys
import pathogenseq as ps
import json

infile = sys.argv[1]
db_prefix = sys.argv[2]
outfile = sys.argv[3]

db_file = db_prefix+".json"
bed_file = db_prefix+".bed"

lt2gene = {}
for l in open(bed_file):
	#Chromosome      5240    7267    Rv0005  gyrB    FLUOROQUINOLONES
	row = l.rstrip().split()
	lt2gene[row[3]] = row[4]

results = json.load(open(infile))
new_results = {"small_variants_dr":[],"small_variants_other":[],"del":[],"lineage":[]}

tmp = ps.db_compare(results,db_file)
dr_muts = [x for x in tmp["variants"] if "annotation" in x]
other_muts = [x for x in tmp["variants"] if "annotation" not in x]
for var in dr_muts:
	#var = {u'genome_pos': 1674782, u'gene_id': u'Rv1484', u'sample': u'por1A', u'chr': u'Chromosome', u'freq': 1.0, u'type': u'missense', 'annotation': [u'ISONIAZID', u'ETHIONAMIDE']
	var["gene"] = lt2gene[var["gene_id"]]
	for drug in var["annotation"]:
		var["drug"] = drug
		new_results["small_variants_dr"].append(json.loads(json.dumps(var)))
for var in other_muts:
	var["gene"] = lt2gene[var["gene_id"]]
	new_results["small_variants_other"].append(var)


json.dump(new_results,open(outfile,"w"))
