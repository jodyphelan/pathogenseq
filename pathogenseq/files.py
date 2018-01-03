import gzip
import os.path
import subprocess

def filecheck(filename):
	"""
	Check if file is there and quit if it isn't
	"""
	if not os.path.isfile(filename):
		print "Can't find %s" % filename
		quit()
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

def run_cmd(cmd,verbose=True):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	cmd = "set -o pipefail; " + cmd
	if verbose==2:
		print "\nRunning command:\n%s" % cmd
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		print "\nRunning command:\n%s" % cmd
		stderr = open("/dev/null","w")
	else:
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr)
	stderr.close()
	if res!=0:
		print "Command Failed! Please Check!"
		quit()

def index_bam(self):
	"""
	Indexing a bam file
	"""
	if os.path.isfile("%(bamfile)s.bai" % self.params):
		cmd = "%(samtools)s index %(bamfile)s" % self.params
		run_cmd(cmd)

def verify_fq(filename):
	"""
	Return True if input is a valid fastQ file
	"""
	FQ = open(filename) if filename[-3:]!=".gz" else gzip.open(filename)
	l1 = FQ.readline()
	if l1[0]!="@":
		print "First character is not \"@\"\nPlease make sure this is fastq format\nExiting..."
		quit()
	else:
		return True

def rm_files(x,verbose=True):
	"""
	Remove a files in a list format
	"""
	for f in x:
		if verbose: print "Removing %s" % f
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
