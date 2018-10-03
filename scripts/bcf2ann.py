import sys
import pathogenseq as ps

bcf = ps.bcf(sys.argv[1])
bcf.get_snp_ann()
