import json
from .files import *
from .mvcf import *
from .fastq import *
from .bam import *
from .mutation_db import *
def profiler(conf_file,prefix,r1=None,r2=None,bam_file=None,call_method="low",min_depth=10,platform="Illumina",mapper="bwa",threads=4):
        conf = json.load(open(conf_file))
        for f in conf:
            filecheck(conf[f])

        if not r1 and not r2 and not bam_file:
            log("Please provide at least one fastQ file (-1) or a BAM file (-b)", True)
        elif (r1 or r2) and bam_file:
            log("Please provide fastQ files or a BAM file but not both",True)
        elif not r1 and r2:
            log("Only second fastQ file provided. If profiling a single ended run, just use '-1'",True)
        if r1 or r2:
            if r2:
                fastq_obj = fastq("bam/"+prefix,conf["ref"],r1,r2,threads=threads)
            else:
                fastq_obj = fastq("bam/"+prefix,conf["ref"],r1,threads=threads)
            if platform=="Illumina":
                bam_obj = fastq_obj.illumina(mapper=mapper)
            elif platform=="minION":
                bam_obj = fastq_obj.minION()
            bam_obj.prefix = "vcf/"+prefix
            bam_obj.params["prefix"] = "vcf/"+prefix
        else:
            log("Using %s.\nPlease ensure that this BAM was made using the same reference as in the database. If unsure what version was used it is best to remap the reads." % bam)
            bam_obj = bam(bam_file,"vcf/"+prefix,conf["ref"])


        bcf = bam_obj.call_variants(call_method=call_method,gff_file=conf["gff"],bed_file=conf["bed"],mixed_as_missing=False,threads=threads,min_dp=min_depth)
        csq = bcf.load_csq_alt(ann_file=conf["ann"],changes=True)
        tmp_bcf = "%s.missing.bcf" % prefix
        missing_pos = get_missing_positions(tmp_bcf)

        results = {"variants":[],"missing":missing_pos}
        for sample in csq:
            results["variants"]  = csq[sample]

        mutations = bam_obj.get_bed_gt(conf["barcode"])
        barcode_mutations = barcode(mutations,conf["barcode"])
        results["barcode"] = barcode_mutations
        results = db_compare(db_file=conf["json_db"],mutations=results)
        return results
