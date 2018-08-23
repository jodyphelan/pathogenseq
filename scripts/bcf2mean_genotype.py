import pathogenseq as ps
import argparse

def main(args):
    bcf = ps.bcf(args.bcf)
    bcf.get_mean_genotype(args.outfile)



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf', help='BCF file')
parser.add_argument('--outfile','-o',default=None,type=str,help='Output file name')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
