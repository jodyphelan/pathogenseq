import json

def db_compare(mutations,db_file):
	db = json.load(open(db_file))
	annotated_mutations = mutations
	for i in range(len(mutations["variants"])):
		#var = {'genome_pos': 6140, 'gene_id': 'Rv0005', 'chr': 'Chromosome', 'freq': 0.975609756097561, 'type': 'missense', 'change': '301V>301L'}
		var = mutations["variants"][i]
		if var["gene_id"] in db and var["change"] in db[var["gene_id"]]:
			annotated_mutations["variants"][i]["annotation"] = db[var["gene_id"]][var["change"]]
	return annotated_mutations
