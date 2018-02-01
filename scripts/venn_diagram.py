#! /usr/bin/env python
import sys
import pathogenseq as ps

bcf_file = sys.argv[1]
bcf = ps.bcf(bcf_file)
bcf.get_venn_diagram_data(sys.argv[2])
