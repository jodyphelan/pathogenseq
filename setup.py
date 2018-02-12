from distutils.core import setup

setup(
	name="pathogenseq",
	version="0.1dev",
	packages=["pathogenseq",],
	license="MIT",
	long_description="Pathogenseq variant calling pipeline",
	scripts=
		['scripts/bam_report.py',
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
		 'scripts/sambamba_depth.py',
		 'scripts/splitchr.py',
		 'scripts/venn_diagram.py',
		 'scripts/fasta2vcf.py'
		]
)
