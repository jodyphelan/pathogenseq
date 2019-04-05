#! /usr/bin/env python
import pathogenseq as ps

ref = sys.argv[1]
query =sys.argv[2]
prefix = sys.argv[3]
ps.mauve_call_variants(ref,query,prefix)
