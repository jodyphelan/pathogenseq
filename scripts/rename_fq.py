import csv
import os
import pathogenseq.files as psf
infile = sys.argv[1]

for row in csv.DictReader(open(infile)):
	f1 = "%s.fastq.gz" % row["Barcode"]
	f2 = "%s.%s.fastq.gz" % (row["Name"],row["Barcode"])
	if psf.filecheck(f1,True)
	os.rename(f1,f2)
