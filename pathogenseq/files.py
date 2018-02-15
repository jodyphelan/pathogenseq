import sys
import gzip
import os.path
import subprocess
import csv
from collections import defaultdict
import json

def init_params():
	conf = json.load(open("%s/%s" % (sys.prefix,"pathogenseq.conf")))
	return conf

def load_tsv(filename):
	meta = {}
	for row in csv.DictReader(open(filename),delimiter="\t"):
		if "sample" not in row:
			print("No sample column...Exiting")
			quit(1)
		meta[row["sample"]] = {}
		columns = set(row)-set(["sample"])
		for c in columns:
			meta[row["sample"]][c] = row[c]
	return columns,meta

def load_bed(filename,columns,key1,key2=None):
	results = defaultdict(lambda: defaultdict(tuple))
	for l in open(filename):
		row = l.rstrip().split()
		if key2:
			if max(columns+[key1,key2])>len(row):
				print "Can't access a column in BED file. The largest column specified is too big"
				quit(1)
			results[row[key1-1]][row[key2-1]] = tuple([row[int(x)-1] for x in columns])
		else:
			if max(columns+[key1])>len(row):
				print "Can't access a column in BED file. The largest column specified is too big"
				quit(1)
			results[row[key1-1]]= tuple([row[int(x)-1] for x in columns])
	return results


def filecheck(filename):
	"""
	Check if file is there and quit if it isn't
	"""
	if not os.path.isfile(filename):
		print("Can't find %s" % filename)
		exit(1)
	else:
		return True

def nofile(filename):
	"""
	Return True if file does not exist
	"""
	if not os.path.isfile(filename):
		return True
	else:
		return False

def bwa_index(ref):
	"""
	Create BWA index for a reference
	"""
	cmd = "bwa index %s" % ref
	run_cmd(cmd)

def run_cmd(cmd,verbose=1):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		print("\nRunning command:\n%s" % cmd)
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		print("\nRunning command:\n%s" % cmd)
		stderr = open("/dev/null","w")
	else:
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr)
	stderr.close()
	if res!=0:
		print("Command Failed! Please Check!")
		exit(1)

def index_bam(bamfile,threads=4):
	"""
	Indexing a bam file
	"""
	if filecheck(bamfile) and nofile(bamfile+".bai"):
		cmd = "samtools index -@ %s %s" % (threads,bamfile)
		run_cmd(cmd)

def index_bcf(bcf,threads=4):
	"""
	Indexing a bcf file
	"""
	if filecheck(bcf) and  nofile(bcf+".csi"):
		cmd = "bcftools index --threads %s %s" % (threads,bcf)
		run_cmd(cmd)
def verify_fq(filename):
	"""
	Return True if input is a valid fastQ file
	"""
	FQ = open(filename) if filename[-3:]!=".gz" else gzip.open(filename)
	l1 = FQ.readline()
	if l1[0]!="@":
		print("First character is not \"@\"\nPlease make sure this is fastq format\nExiting...")
		exit(1)
	else:
		return True

def rm_files(x,verbose=True):
	"""
	Remove a files in a list format
	"""
	for f in x:
		if verbose: print("Removing %s" % f)
		os.remove(f)

def file_len(filename):
	"""
	Return length of a file
	"""
	filecheck(filename)
	for l in subprocess.Popen("wc -l %s" % filename,shell=True,stdout=subprocess.PIPE).stdout:
		res = l.rstrip().split()[0]
	return int(res)

def gz_file_len(filename):
	"""
	Return lengths of a gzipped file
	"""
	filecheck(filename)
	for l in subprocess.Popen("gunzip -c %s |wc -l" % filename,shell=True,stdout=subprocess.PIPE).stdout:
		res = l.rstrip().split()[0]
	return int(res)

def download_from_ena(acc):
	if len(acc)==9:
		dir1 = acc[:6]
		cmd = "wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s*" % (dir1,acc,acc)
	elif len(acc)==10:
		dir1 = acc[:6]
		dir2 = "00"+acc[-1]
		cmd = "wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s*" % (dir1,dir2,acc,acc)
	else:
		print("Check Accession: %s" % acc)
		exit(1)
	run_cmd(cmd)

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def programs_check(programs):
	for p in programs:
		if which(p)==None:
			print "Can't find %s in path... Exiting." % p
			quit(1)
