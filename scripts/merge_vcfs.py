#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse


def main(args):
	vcf = ps.vcf_merge(args.samples,args.ref,args.prefix,args.mappability_filter,args.mappability_file,args.vcf_dir,args.min_dp,args.keep_samples,args.fmiss,args.miss_cut,args.mix_cut,args.low_cov,args.bed_include,args.bed_exclude,args.threads,args.vcf_ext)
	vcf.merge()
	vcf.extract_variants()
	vcf.filt_non_uniq()
	vcf.sample_filt()
	vcf.mask_mixed()
	bcf = vcf.get_bcf_obj()
	bcf.vcf_to_fasta(args.prefix+".snps.fa",args.threads)



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('samples', help='First read file')
parser.add_argument('ref', help='Second read file')
parser.add_argument('prefix', help='Reference Sequence')
parser.add_argument('--threads','-t', type=int, default=1, help='Number of threads')
parser.add_argument('--mappability_filter',action="store_true", help='Reference Sequence')
parser.add_argument('--mappability_file','-m', help='Reference Sequence')
parser.add_argument('--min_dp', default=10, type=int, help='Reference Sequence')
parser.add_argument('--keep_samples', default=None,help='Reference Sequence')
parser.add_argument('--fmiss', default=0.1,type=float,help='The maximum fraction of missing data to keep SNP position')
parser.add_argument('--miss_cut', default=0.15,type=float,help='Reference Sequence')
parser.add_argument('--mix_cut',default=0.15,type=float, help='Reference Sequence')
parser.add_argument('--low_cov', default=False,help='Reference Sequence')
parser.add_argument('--bed_include',default=None, help='Reference Sequence')
parser.add_argument('--bed_exclude', default=None,help='Reference Sequence')
parser.add_argument('--vcf_dir', default=".",help='Reference Sequence')
parser.add_argument('--vcf_ext', default="vcf.gz",help='Reference Sequence')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
