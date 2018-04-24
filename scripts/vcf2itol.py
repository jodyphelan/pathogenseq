#! /usr/bin/env python
import sys
import pathogenseq as ps

x = ps.bcf(sys.argv[1])
x.itol_from_bcf(sys.argv[2])
