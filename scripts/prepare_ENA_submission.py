import sys
import csv
import argparse
import subprocess

def main(args):
	sample_csv = "%s.samples.tsv" % args.prefix
	run_csv = "%s.runs.tsv" % args.prefix
	sample_fieldnames = ["sample_alias","tax_id","scientific_name","common_name","sample_title","sample_description"]
	SAMP = open(sample_csv,"w")
	SAMP.write("#checklist_accession\tERC000011\n#unique_name_prefix\n")

	samp_writer = csv.DictWriter(SAMP,fieldnames=sample_fieldnames,delimiter="\t")
	samp_writer.writeheader()

	RUNS = open(run_csv,"w")
	runs_fieldnames = ["project_accession", "project_alias", "sample_alias", "experiment_alias", "run_alias", "library_name", "library_source", "library_selection", "library_strategy", "design_description", "library_construction_protocol", "instrument_model", "file_type", "library_layout", "insert_size", "forward_file_name", "forward_file_md5", "forward_file_unencrypted_md5", "reverse_file_name", "reverse_file_md5", "reverse_file_unencrypted_md5"]
	run_writer = csv.DictWriter(RUNS,fieldnames=runs_fieldnames,delimiter="\t")
	run_writer.writeheader()
	md5 = {}

	subprocess.call("awk -F ',' '$1!=\"id\" {print \"%s/\"$2\",%s\"$3}' %s | tr ',' '\\n' |  parallel -j 4 \" if [ ! -f {}.md5 ]; then md5sum {} > {}.md5; fi;\"" % (args.fastq_dir,args.fastq_dir,args.sample_file), shell=True)
	print("awk -F ',' '$1!=\"id\" {print \"%s/\"$2\",%s\"$3}' %s | tr ',' '\\n' |  parallel -j 4 \" if [ ! -f {}.md5 ]; then md5sum {} > {}.md5; fi;\"")

	for row in csv.DictReader(open(args.sample_file)):
		print(row["id"])
		md5[row["R1"]] = open("%s/%s.md5" % (args.fastq_dir,row["R1"])).readline().strip().split()[0]
		md5[row["R2"]] = open("%s/%s.md5" % (args.fastq_dir,row["R2"])).readline().strip().split()[0]
		tmp_samp = {"sample_alias":args.sample_alias_prefix+"_"+row["id"],"tax_id":args.tax_id,"scientific_name":args.scientific_name,"common_name":"","sample_title":args.sample_title,"sample_description":""}
		samp_writer.writerow(tmp_samp)
		tmp_run = {"project_accession":args.project, "project_alias":"", "sample_alias":args.sample_alias_prefix+"_"+row["id"], "experiment_alias":"", "run_alias":"", "library_name":"unspecified", "library_source":args.library_source, "library_selection":args.library_selection, "library_strategy":args.library_strategy, "design_description":"", "library_construction_protocol":"", "instrument_model":args.instrument_model, "file_type":"fastq", "library_layout":"PAIRED", "insert_size":args.insert_size, "forward_file_name":row["R1"], "forward_file_md5":md5[row["R1"]], "forward_file_unencrypted_md5":"", "reverse_file_name":row["R2"], "reverse_file_md5":md5[row["R2"]], "reverse_file_unencrypted_md5":""}
		run_writer.writerow(tmp_run)

	SAMP.close()
	RUNS.close()
parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('prefix', help='file prefix')
parser.add_argument('sample_file', help='BAM file')
parser.add_argument('sample_title', help='Reference Sequence')
parser.add_argument('scientific_name', help='BAM file')
parser.add_argument('tax_id', help='Prefix for files')
parser.add_argument('project', help='Prefix for files')
parser.add_argument('--library_source',default="GENOMIC",choices=["GENOMIC","TRANSCRIPTOMIC","METAGENOMIC"], help='Prefix for files')
parser.add_argument('--library_selection',default="PCR")
parser.add_argument('--library_strategy', choices=["WGS","WGA","RNA-seq"],default="WGS")
parser.add_argument('--instrument_model',default="Illumina HiSeq 2500", choices=["Illumina HiSeq 2000","Illumina HiSeq 2500","Illumina HiSeq 3000","Illumina HiSeq 4000","Illumina MiSeq","Illumina NovaSeq 6000"])
parser.add_argument('--insert_size',default="500")
parser.add_argument('--sample_alias_prefix',default="")
parser.add_argument('--fastq_dir',default=".")
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
