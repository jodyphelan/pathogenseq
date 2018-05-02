#! /usr/bin/env python
import sys
import pathogenseq as ps

bcf_file,ref_file = sys.argv[1:]
bcf = ps.bcf(bcf_file)
bcf.generate_consensus(ref_file)
