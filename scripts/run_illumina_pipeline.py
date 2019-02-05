#! /usr/bin/env python
import sys
import pathogenseq as ps
import csv
import argparse
import json
import statistics

def process(args):
	out_script = "%s.run.sh" % args.prefix
	O = open(out_script,"w")
	ps.filecheck(args.ref)

	samples = []
	args.sample_file = "%s.samples.txt" % args.prefix
	for row in csv.DictReader(open(args.csv)):
		args.id = row["ID"]
		samples.append(row["ID"])
		args.r1 = "%s/%s" % (args.fastq_dir,row["R1"])
		ps.filecheck(args.r1)
		#params["centrifuge"] = "--centrifuge %s" % args.centrifuge if args.centrifuge else ""
		args.r2 = "%s/%s" % (args.fastq_dir,row["R2"])
		ps.filecheck(args.r2)
		O.write("illumina_pipeline.py %(ref)s %(r1)s %(r2)s %(id)s -t %(threads)s -m %(mapper)s \n" % vars(args))
	O.close()
	open(args.sample_file,"w").write("\n".join(samples))
	if not args.dry:
		ps.run_cmd("bash %s" % out_script)


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('process', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('csv', help='First read file')
parser_sub.add_argument('ref', help='First read file')
parser_sub.add_argument('prefix', help='First read file')
parser_sub.add_argument('--platform','-m',choices=["illumina","minION"],default="illumina", help='First read file')
parser_sub.add_argument('--fastq_dir','-f',default=".",type=str, help='First read file')
parser_sub.add_argument('--threads',"-t",type=int,default=1, help='First read file')
parser_sub.add_argument('--mapper',type=str,choices=["bwa","minimap2","bowtie2"],default="bwa", help='First read file')
#parser_sub.add_argument('--centrifuge','-c',type=str,default=None)
parser_sub.add_argument('--window',default=None,type=int, help='First read file')
parser_sub.add_argument('--dry',action="store_true", help='First read file')
parser_sub.set_defaults(func=process)

parser_sub = subparsers.add_parser('collect', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('prefix', help='First read file')
parser_sub.add_argument('--samples','-s', help='First read file')
parser_sub.add_argument('--dist',default=10,type=int, help='First read file')
parser_sub.set_defaults(func=collect)


args = parser.parse_args()
args.func(args)
