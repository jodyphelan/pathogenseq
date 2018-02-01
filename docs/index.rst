.. pathogenseq documentation master file, created by
   sphinx-quickstart on Thu Dec 21 16:11:15 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pathogenseq pipeline documentation
=======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
This package contains the code to run many of the basic QC and processing operations on NGS data.

Installation
------------

.. code-block:: bash

    git clone git@github.com:jodyphelan/pathogenseq.git
    cd  pathogenseq
    python setup.py install

Package structure
-----------------
The package is divided into modules:

1. *qc* - for performing QC operations
2. *map_call_snps* - performing mapping and calling of SNPs
3. *vcf_merge* - merging gVCF into multi-sample VCFs
4. *mvcf* - common operations on multi-sample VCFs
5. *fasta* - representing a fasta as a python data structure
6. *files* - functions to handle calls to the commandline and ckecking files are present
7. *logger* - **(in development)** will be used to create log files

These modules are imported into the main *pathogenseq* module and their associated classes and functions can be accessed like this:

>>> import pathogenseq as ps
>>> fastqqc = ps.qc_fastq("prefix","/path/to/read_1","/path/to/read_2")

The idea is that the modules are highly flexible and are easy to integrate into scripts.
Some example scripts are present in the ``scripts`` directory.

A very simple example script to calculate fraction of genome covered at predefined coverage cutoff could be the following:

.. code-block:: python

		import sys
		import pathogenseq as ps

		bam_file = sys.argv[1]
		ref_file = sys.argv[2]

		bamqc = ps.qc_bam(bam_file,ref_file,[1,5,10])
		print "Fraction genome coverage at depth:"
		print "1 = %s" % (bamqc.genome_cov[1])
		print "5 = %s" % (bamqc.genome_cov[5])
		print "10 = %s" % (bamqc.genome_cov[10])




Performing QC on a fastQ file(s)
--------------------------------

To perform QC on a fastq files we can use the ``qc_fastq`` class. To initialise the class we need the location of the read files and a prefix for any output files:

>>> import pathogenseq as ps
>>> fastqqc = ps.qc_fastq("prefix","/path/to/read_1.fq.gz","/path/to/read_2.fq.gz")

We can now calculate metrics such as the number of reads and median read length.

>>> # number of reads
>>> fastqqc.read_num
>>> # median read len
>>> fastqqc.median_read_len
>>> # get approx depth for a 4.4 Mb genome
>>> fastqqc.approx_depth("4.4Mb")

Additionally we can run kraken (using a pre-generated database). We can also create filtered fastq files with only reads corresponding to our genome of interest (using the NCBI Taxonomic ID: https://www.ncbi.nlm.nih.gov/taxonomy). Multiple taxonomies can be supplied as a comma delimited string.

>>> # filter out any reads not belonging to the Mycobacterium tuberculosis complex
>>> # this will create two files with the user provided prefix and a .kraken_filt.fastq.gz extension
>>> fastqqc.run_kraken("/opt/storage2/ernest/kraken/kraken/standard_db","77643,1773,78331,33894,1765")
>>> # get a mapping class object
>>> fastqqc.get_mapper_from_kraken("/path/to/reference")

Performing mapping
------------------
We can use the ``mapping`` class to access basic methods involved with QC, mapping and variant detections.

>>> import pathogenseq as ps
>>> mapper = ps.fastq("/path/to/read_1.fq.gz","/path/to/read_2.fq.gz","/path/to/reference.fasta","prefix",threads=4)
>>> # perform basic trimming
>>> mapper.trim()
>>> # perform mapping with BWA (bowtie available too)
>>> mapper.map()
>>> # get qc_bam class from output bam
>>> mapper.get_bam_qc()

Extract stats on BAM files
-----------------------------

We can use the ``qc_bam`` class to access methods to extract stats from bam files.

>>> import pathogenseq as ps
>>> bamqc = ps.qc_bam("/path/to/bam","/path/to/ref.fasta")
>>> # get a dict with the fraction of genome covered for predefined depth thresholds
>>> bamqc.genome_cov
>>> # get median depth
>>> bamqc.med_dp
>>> # get percentage of reads mapping
>>> bamqc.pct_reads_mapped
>>> # get a dict with the coverage across each position in the genome
>>> # keys are the difference chromosomes
>>> bamqc.ref_dp
>>> # dump the ref_dp dict into a json file
>>> bamqc.save_cov("/path/to/output.json")
>>> # get the mean coverage across regions in a BED file
>>> bamqc.bed_cov("/path/to/bed")
>>> # create coverage plot for a chromosome
>>> bamqc.plot_cov("Chromosome","/path/to/output.png")

Merging gVCF files
------------------

Merging vcf files can be done with the ``vcf_merge`` class. VCF files must be **bgzipped**
If you have a file containing the sample names (one per line) and directory structure like this:

.. code-block:: bash

    └── vcf
        ├── sample1.vcf.gz
        ├── sample2.vcf.gz
        └── sample3.vcf.gz

>>> import pathogenseq as ps
>>> vcf = vcf_merge("/path/to/samples.txt","/path/to/ref.fasta","prefix")
>>> # perform merging
>>> vcf.merge()
>>> # extract only high quality SNPs from merged file
>>> vcf.extract_variants()
>>> # filter out non unique regions of the genome
>>> vcf.filt_non_uniq()
>>> # filter out samples
>>> vcf.sample_filt()
>>> # mask mixed positions
>>> vcf.mask_mixed()
>>> # generate a multi-fasta file with variants
>>> vcf.generate_fasta()

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
