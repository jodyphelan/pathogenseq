import pathogenseq
import argparse









def main(args):
	x= pathogenseq.mapping(args.r1,args.r2,args.ref,args.prefix,threads=args.threads)
	x.trim()
	x.map()
	x.call_snps()






parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--r1','-1', help='First read file')
parser.add_argument('--r2','-2', help='Second read file')
parser.add_argument('--ref','-r', help='Reference Sequence')
parser.add_argument('--threads','-t', type=int, default=1, help='Number of threads')
parser.add_argument('--prefix','-p', help='Prefix for files')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
