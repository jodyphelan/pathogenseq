#! /usr/bin/env python
import sys
import pathogenseq as ps
import argparse


def main(args):
	if args.mappability_file:
		args.mappability_filter=True
	bcf_variant_pos = "%s.prefilt.bcf" % args.prefix
	bcf_sample_filt = "%s.sample_filt.bcf" % args.prefix
	bcf_uniq_filt = "%s.uniq_filt.bcf" % args.prefix
	bcf_variant_filt = "%s.variant_filt.bcf" % args.prefix
	bcf_masked_filt = "%s.mix_masked.bcf" % args.prefix
	fasta_snps = "%s.snps.fa" % args.prefix

	merged = ps.vcf_merge(args.samples,args.ref,args.prefix,args.vcf_dir,args.vcf_ext,args.threads,args.min_dp).merge()

	vcf = merged.extract_variants(bcf_variant_pos,min_dp=args.min_dp,bed_include=args.bed_include,bed_exclude=args.bed_exclude)

	if args.mappability_filter:
		if not args.mappability_file:
			create_mappability_file(args.ref,args.threads)
			args.mappability_file = "%s.genome.mappability.bed" % args.prefix
		vcf = vcf.filt_non_uniq(args.mappability_file,bcf_uniq_filt)

	vcf = vcf.sample_filt(bcf_sample_filt,miss_cut=args.miss_cut,mix_cut=args.mix_cut,keep_samples=args.keep_samples)

	vcf = vcf.filt_variants(bcf_variant_filt,fmiss=args.fmiss,threads=args.threads)

	vcf = vcf.mask_mixed(bcf_masked_filt)

	vcf.vcf_to_fasta_alt(outfile=fasta_snps,ref_file=args.ref,threads=args.threads)





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
parser.add_argument('--bed_include',default=None, help='Reference Sequence')
parser.add_argument('--bed_exclude', default=None,help='Reference Sequence')
parser.add_argument('--vcf_dir', default=".",help='Reference Sequence')
parser.add_argument('--vcf_ext', default="gbcf",help='Reference Sequence')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
