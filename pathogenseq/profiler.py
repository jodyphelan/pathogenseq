from __future__ import division
import json
from .files import *
from .mvcf import *
from .fastq import *
from .bam import *
from .mutation_db import *
from .qc import *
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
                fastq_obj = fastq(prefix,conf["ref"],r1,r2,threads=threads)
            else:
                fastq_obj = fastq(prefix,conf["ref"],r1,threads=threads)
            if platform=="Illumina":
                bam_obj = fastq_obj.illumina(mapper=mapper)
            elif platform=="minION":
                bam_obj = fastq_obj.minION()
        else:
            log("Using %s\nPlease ensure that this BAM was made using the same reference as in the database.\nIf you are not sure what reference was used it is best to remap the reads." % bam_file)
            bam_obj = bam(bam_file,prefix,conf["ref"],platform=platform)
        bcf_obj = bam_obj.call_variants(call_method=call_method,gff_file=conf["gff"],bed_file=conf["bed"],mixed_as_missing=False if platform == "Illumina" else True,threads=threads,min_dp=min_depth)
        #bcf_obj = bcf(prefix+".csq.bcf")
        csq = bcf_obj.load_csq_alt(ann_file=conf["ann"],changes=True)
        bamqc = qc_bam(bam_obj.bam_file,conf["ref"],bed_file=conf["bed"])

        tmp_bcf = "%s.missing.bcf" % prefix
        missing_pos = get_missing_positions(tmp_bcf)
        bed_pos = load_bed(conf["bed"],[1,2,3,4],4,intasint=True)
        miss_gene = {}
        for gene in bed_pos:
            miss_pos = 0
            chrom,start,end = bed_pos[gene][0:3]
            start = int(start)
            end = int(end)
            for i in range(start,end):
                if [chrom,i] in missing_pos or bamqc.ref_dp[chrom][i]==0:
                     miss_pos+=1
            miss_gene[gene] = miss_pos/(end-start)
        results = {"variants":[],"missing":miss_gene,"qc":{"pct_reads_mapped":bamqc.pct_reads_mapped,"num_reads_mapped":bamqc.num_reads_mapped}}
        for sample in csq:
            results["variants"]  = csq[sample]

        mutations = bam_obj.get_bed_gt(conf["barcode"])
        barcode_mutations = barcode(mutations,conf["barcode"])
        results["barcode"] = barcode_mutations
        results = db_compare(db_file=conf["json_db"],mutations=results)
        return results
