#! /usr/bin/env python
import sys
import pathogenseq.files as ps
import csv
import argparse

def main(args):
	out_script = "%s.run.sh" % args.prefix
	O = open(out_script,"w")
	samples = []
	for row in csv.DictReader(open(args.sample_file)):
		params = {}
		params["ref_file"] = "%s/%s" % (args.ref_dir,row["Reference"])
		ps.filecheck(params["ref_file"])
		params["r1"] = "%s/%s" % (args.fastq_dir,row["ReadF"])
		ps.filecheck(params["r1"])
		params["prefix"] = row["ID"]
		samples.append(row["ID"])
		params["threads"] = args.threads
		params["mapper"] = args.mapper
		params["centrifuge"] = "--centrifuge %s" % args.centrifuge if args.centrifuge else ""

		if "Primers" in row and row["Primers"]!="NA":
			params["primers"] = "--primers %s/%s" % (args.primer_dir,row["Primers"])
		else:
			params["primers"] = ""
		if args.platform=="illumina":
			params["r2"] = "%s/%s" % (args.fastq_dir,row["ReadR"])
			ps.filecheck(params["r2"])
			O.write("illumina_pipeline.py %(ref_file)s %(r1)s %(r2)s %(prefix)s -t %(threads)s -m %(mapper)s %(primers)s %(centrifuge)s\n" % params)
		else:
			params["window"] = "--window %s" % args.window if args.window else ""
			O.write("minION_pipeline.py %(ref_file)s %(r1)s %(prefix)s -t %(threads)s %(primers)s %(centrifuge)s %(window)s\n" % params)
	O.close()
	open("%s.samples.txt" % args.prefix,"w").write("\n".join(samples))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('sample_file', help='First read file')
parser.add_argument('prefix', help='First read file')
parser.add_argument('--platform','-m',choices=["illumina","minION"],default="illumina", help='First read file')
parser.add_argument('--ref_dir','-r',default=".",type=str, help='First read file')
parser.add_argument('--fastq_dir','-f',default=".",type=str, help='First read file')
parser.add_argument('--primer_dir',"-p",default=".",type=str, help='First read file')
parser.add_argument('--threads',"-t",type=int,default=1, help='First read file')
parser.add_argument('--mapper',type=str,choices=["bwa","minimap2","bowtie2"],default="bwa", help='First read file')
parser.add_argument('--centrifuge','-c',type=str,default=None)
parser.add_argument('--window',default=None,type=int, help='First read file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
