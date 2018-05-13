import json

def db_compare(mutations,db_file):
	db = json.load(open(db_file))
	for var in mutations:
		#var = {'genome_pos': 6140, 'gene_id': 'Rv0005', 'chr': 'Chromosome', 'freq': 0.975609756097561, 'type': 'missense', 'change': '301V>301L'}
		if var["gene_id"] in db and var["change"] in db[var["gene_id"]]:
			print "%s %s %s" % (var["sample"],var["gene_id"],var["change"])
