from distutils.core import setup, Command
import sys
import os
import subprocess
import json

default_conf = {
	"raxmlHPC":"raxmlHPC",
	"examl":"examl",
	"parse_examl":"parse-examl",
	"bwa":"bwa",
	"bowtie2":"bowtie2",
	"bowtie2-build":"bowtie2-build",
	"minimap2":"minimap2",
	"tabix":"tabix",
	"bgzip":"bgzip",
	"bcftools":"bcftools",
	"samtools":"samtools",
	"htsbox":"htsbox"
}

class InstallCommand(Command):
    description = "Installs the foo."
    user_options = [
        ('arch=', "a", 'Specify the arch'),
    ]
    def initialize_options(self):
        self.arch = "Linux"
    def finalize_options(self):
        assert self.arch in ('Linux', 'osX'), 'Invalid Arch!'
    def run(self):
		subprocess.call("bash install_prerequisites.sh",shell=True)
		cwd = os.path.dirname(os.path.realpath(__file__))
		new_conf = default_conf
		for d in new_conf:
			new_conf[d] = "%s/bin/%s" % (cwd,default_conf[d])
		json.dump(new_conf,open("pathogenseq.conf","w"))

class init_conf(Command):
    description = "Installs the foo."
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
		json.dump(default_conf,open("pathogenseq.conf","w"))

setup(

	name="pathogenseq",
	cmdclass={
		"install_prerequisites": InstallCommand,
		"init_conf": init_conf
	},
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
		 'scripts/rename_fq.py'
		],
	data_files = [(sys.prefix,["pathogenseq.conf"])],

)
