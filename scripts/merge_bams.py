import pathogenseq as ps
import sys
import argparse

def merge_sample(sample,cram=False,reference=None,threads=20):
	cram_flag = "-CT %s" % reference if cram else ""
	extension = ".cram" if cram else ".bam"
	samples = sample.split("_")

	for x in samples:
		ps.filecheck(x+extension)
	bams = " ".join([x+extension for x in samples])
	first_bam = samples[0]+extension
	temp_bam = sample+".tmp.bam"
	final_bam = sample+extension
	header = sample+".header"
	cmd = "samtools view %s -H | sed 's/%s/%s/g' > %s" % (first_bam,samples[0],sample,header)
	ps.run_cmd(cmd)
	cmd = "samtools merge -@ %s %s %s" % (threads,temp_bam,bams)
	ps.run_cmd(cmd)
	cmd = "samtools reheader %s %s | samtools view -@ %s %s > %s" % (header,temp_bam,threads,cram_flag,final_bam)
	ps.run_cmd(cmd)
	cmd = "rm %s %s" % (header,temp_bam)
	ps.run_cmd(cmd)

def main(args):
	if not args.samples_file and not args.sample:
		ps.log("Provide either --sample or --samples_file... Exiting",ext=True)
	if args.samples_file and args.sample:
		ps.log("Provide either --sample or --samples_file but not both... Exiting",ext=True)
	if args.cram and not args.reference:
		ps.log("Provide reference for cram compression... Exiting",ext=True)
	if args.sample:
		merge_sample(args.sample,args.cram,args.reference)
	else:
		for l in open(args.samples_file):
			merge_sample(l.rstrip(),args.cram,args.reference)





parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples_file','-f',default=None,help='Sample file')
parser.add_argument('--sample','-s',default=None,help='Individual sample')
parser.add_argument('--reference','-r',type=str,default=None,help='Cram file')
parser.add_argument('--cram','-c',action="store_true",help='Cram file')
parser.add_argument('--threads','-t',default=20,type=int,help='Cram file')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
