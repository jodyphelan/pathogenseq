#! /usr/bin/env python
import sys
import pathogenseq.files as ps
import csv
import argparse

def main(args):
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
		O.write("tb-profiler profile -p %(id)s -a %(id)s.bam -t %(threads)s\n" % vars(args))

	O.write("merge_vcfs.py %(sample_file)s %(ref)s %(prefix)s\n" % vars(args))
	args.snps_aln_file = "%s.snps.fa" % args.prefix
	O.write("raxml-ng --msa %(snps_aln_file)s --search --model GTR+G --threads `raxml-ng --msa %(snps_aln_file)s --parse --model GTR+G | grep MPI | awk '{print $9}'`\n" % vars(args))
	O.write("tb-profile collate\n")



	O.close()
	open(args.sample_file,"w").write("\n".join(samples))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('csv', help='First read file')
parser.add_argument('ref', help='First read file')
parser.add_argument('prefix', help='First read file')
parser.add_argument('--platform','-m',choices=["illumina","minION"],default="illumina", help='First read file')
parser.add_argument('--fastq_dir','-f',default=".",type=str, help='First read file')
parser.add_argument('--threads',"-t",type=int,default=1, help='First read file')
parser.add_argument('--mapper',type=str,choices=["bwa","minimap2","bowtie2"],default="bwa", help='First read file')
#parser.add_argument('--centrifuge','-c',type=str,default=None)
parser.add_argument('--window',default=None,type=int, help='First read file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
