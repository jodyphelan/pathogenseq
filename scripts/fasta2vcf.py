#! /usr/bin/env python
import sys
import pathogenseq as ps

ref_file = sys.argv[1]
query_file = sys.argv[2]
prefix = sys.argv[3]
ps.mauve_call_variants(ref_file,query_file,prefix)
cmd = "bgzip -f %s.vcf" % prefix
ps.run_cmd(cmd)
