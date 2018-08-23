import pathogenseq as ps
import argparse
import csv

def main(args):
	bcf = ps.bcf(args.bcf)
	#bcf.get_mean_genotype()
	geno_file = bcf.prefix+".geno"
	meta = {}
	for s in bcf.samples:
		meta[s] = {}
	for row in csv.DictReader(open(args.csv)):
		for pheno in row.keys():
			if pheno=="id": continue
			if row['id'] not in meta: continue
			meta[row["id"]][pheno] = row[pheno]
	phenos = [x.rstrip() for x in open(args.phenos).readlines()]
	cmd_file = ps.get_random_file()
	X = open(cmd_file,"w")
	for pheno in phenos:
		pheno_file = "%s.pheno" % pheno
		if pheno not in row:
			ps.log("%s not in CSV file"%pheno,True)
		P = open(pheno_file,"w")
		P.write("\n".join([meta[s][pheno] if pheno in meta[s] else "NA" for s in bcf.samples]))
		P.close()
		X.write("gemma -p %s -g %s -gk 1 -o %s && gemma  -lmm 1 -p %s -g %s  -k output/%s.cXX.txt  -o %s\n" % (pheno_file,geno_file,pheno,pheno_file,geno_file,pheno,pheno))
	X.close()

	if args.preprocess:
		ps.log("Preprocessing finished\n", True)
	else:
		ps.run_cmd("cat %s | parallel -j %s" % (cmd_file,args.threads))
parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf', help='BCF file')
parser.add_argument('csv', help='CSV file')
parser.add_argument('phenos', help='Columns file')
parser.add_argument('--preprocess',default=False,action='store_true', help='Columns file')
parser.add_argument('--threads','-t',default=4,type=int, help='Columns file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
