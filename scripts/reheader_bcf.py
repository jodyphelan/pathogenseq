#! /usr/bin/env python
import pathogenseq as ps
import sys

bcf = ps.bcf(sys.argv[1])
bcf.reheader(sys.argv[2])

