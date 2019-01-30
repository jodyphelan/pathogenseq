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
		O.write("tb-profiler profile -p %(id)s -a %(id)s.bam -t %(threads)s\n" % vars(args))

	O.write("merge_vcfs.py %(sample_file)s %(ref)s %(prefix)s\n" % vars(args))
	args.snps_aln_file = "%s.snps.fa" % args.prefix
	O.write("raxml-ng --msa %(snps_aln_file)s --search --model GTR+G --threads `raxml-ng --msa %(snps_aln_file)s --parse --model GTR+G | grep MPI | awk '{print $9}'`\n" % vars(args))
	O.write("tb-profiler collate\n")
	O.close()
	open(args.sample_file,"w").write("\n".join(samples))
	if not args.dry:
		ps.run_cmd("bash %s" % out_script)

def collect(args):
	samples = [x.rstrip() for x in open("%s.samples.txt" % args.prefix).readlines()]
	dataset_stats = {}
	coverage = []
	for s in samples:
		tmp = json.load(open("%s.stats.json" % s))
		coverage.append(tmp["med_dp"])
	dataset_stats["median_dp"] = statistics.median(coverage)
	dataset_stats["max_dp"] = max(coverage)
	dataset_stats["min_dp"] = min(coverage)
	bcf = ps.bcf("%s.mix_masked.bcf" % args.prefix)
	bcf_stats = bcf.load_stats()
	dataset_stats["num_variants"] = bcf_stats["SN"]["number of records"]
	dists = bcf.get_plink_dist()
	min_dist = 1000000
	for i in range(len(dists)):
		for j in range(len(dists)):
			if j<=i:continue
			if dists[i][j]<min_dist:min_dist = dists[i][j]
	dataset_stats["min_dist"] = min_dist
	dataset_stats["num_snps_one_sample"] = bcf_stats["AF"]["SNP"][sorted(list(bcf_stats["AF"]["SNP"].keys()))[0]]
	dataset_stats["num_snps_one_sample_pct"] = dataset_stats["num_snps_one_sample"]/dataset_stats["num_variants"]
	tbprofiler_json = json.load(open("%s.json" % args.prefix))
	dataset_stats["Drug-resistant"] = 0
	dataset_stats["MDR"] = 0
	dataset_stats["XDR"] = 0
	dataset_stats["Sensitive"] = 0
	for s in tbprofiler_json:
		dataset_stats[tbprofiler_json[s]["drtype"]]+=1



	O = open("%s.stats.txt"%args.prefix,"w")
	for key in dataset_stats:
		O.write("%s\t%s\n" % (key,dataset_stats[key]))
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
parser_sub.set_defaults(func=collect)


args = parser.parse_args()
args.func(args)
