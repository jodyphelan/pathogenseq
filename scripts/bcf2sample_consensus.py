import pathogenseq as ps
import sys

bcf = ps.bcf(sys.argv[1])
bcf.per_sample_bcf2fa(sys.argv[2],sys.argv[3])
