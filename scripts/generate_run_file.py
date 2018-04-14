import sys
import pathogenseq.files as ps
import csv
import argparse

def main(args):

	O = open(args.out_script,"w")
	for row in csv.DictReader(open(args.sample_file)):
		params = {}
		params["ref_file"] = "%s/%s" % (args.ref_dir,row["Reference"])
		ps.filecheck(params["ref_file"])
		params["r1"] = "%s/%s" % (args.fastq_dir,row["ReadF"])
		ps.filecheck(params["r1"])
		params["prefix"] = row["ID"]
		params["threads"] = args.threads
		if "Primers" in row and row["Primers"]!="NA":
			params["primers"] = "--primers %s/%s" % (args.primer_dir,row["Primers"])
		else:
			params["primers"] = ""
		if args.platform=="illumina":
			params["r2"] = "%s/%s" % (args.fastq_dir,row["ReadR"])
			ps.filecheck(params["r2"])
			O.write("illumina_pipeline.py %(ref_file)s %(r1)s %(r2)s %(prefix)s -t %(threads)s %(primers)s\n" % params)
		else:
			O.write("minION_pipeline.py %(ref_file)s %(r1)s %(prefix)s -t %(threads)s %(primers)s\n" % params)
	O.close()

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('sample_file', help='First read file')
parser.add_argument('out_script', help='First read file')
parser.add_argument('--platform','-m',choices=["illumina","minION"],default="illumina", help='First read file')
parser.add_argument('--ref_dir','-r',default=".",type=str, help='First read file')
parser.add_argument('--fastq_dir','-f',default=".",type=str, help='First read file')
parser.add_argument('--primer_dir',"-p",default=".",type=str, help='First read file')
parser.add_argument('--threads',"-t",type=int,default=1, help='First read file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
