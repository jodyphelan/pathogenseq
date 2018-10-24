import pathogenseq as ps
import argparse

def reheader(bamfile,oldname,newname,threads):
	tmpfile = "%s.header" % bamfile
	cmd = "samtools view %s -H | sed 's/%s/%s/g' > %s" % (bamfile,oldname,newname,tmpfile)
	ps.run_cmd(cmd)
	newbamfile = "%s.reheader.bam" % (bamfile.replace(".bam",""))
	cmd = "samtools reheader %s %s | samtools view -@ %s -b > %s" % (tmpfile,bamfile,newbamfile,threads)
	ps.run_cmd(cmd)

def main(args):
	reheader(args.bam,args.oldname,args.newname,args.threads)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bam', help='BCF file')
parser.add_argument('oldname', help='CSV file')
parser.add_argument('newname', help='Columns file')
parser.add_argument('--threads','-t',default=4,type=int, help='Columns file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
