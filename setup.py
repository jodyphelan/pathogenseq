import setuptools
import sys
import os
import subprocess
import json


setuptools.setup(

	name="pathogenseq",
	version="0.1dev",
	packages=["pathogenseq",],
	license="MIT",
	long_description="Pathogenseq variant calling pipeline",
	scripts=
		['scripts/bcf2consensus.py',
		 'scripts/bam_report.py',
		 'scripts/bcf2fasta.py',
		 'scripts/bam2vcf.py',
		 'scripts/collate_json.py',
		 'scripts/cov_plt.py',
		 'scripts/examl.py',
		 'scripts/fastq_report.py',
		 'scripts/mapping.py',
		 'scripts/merge_vcfs.py',
		 'scripts/minION_pipeline.py',
		 'scripts/Mtb_pipeline.py',
		 'scripts/Mtb_pipeline_no_filt.py',
		 'scripts/sambamba_depth.py',
		 'scripts/splitchr.py',
		 'scripts/venn_diagram.py',
		 'scripts/fasta2vcf.py',
		 'scripts/pathogen-profiler.py',
		 'scripts/combine_dict_list.py',
		 'scripts/illumina_pipeline.py',
		 'scripts/generate_run_file.py',
		 'scripts/vcf2itol.py',
		 'scripts/csv2fasta.py',
		 'scripts/rename_fq.py',
		 'scripts/profiler2tbprofiler.py',
		 'scripts/select_reference.py',
		 'scripts/bcf2dist.py',
		 'scripts/bcf_sample_stats.py',
		 'scripts/extract_dosage.py'
		],
)
