#! /usr/bin/env python
import pathogenseq as ps

bcf = ps.bcf(sys.argv[1])
bcf.extract_dosage(sys.argv[2])
