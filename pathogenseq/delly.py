from .files import *
from .mvcf import bcf
class delly_bcf(bcf):
	def __init__(self,filename):
		 bcf.__init__(self,filename)
	def get_robust_calls(self):
		results = []
		for l in cmd_out(" bcftools query -f '%%CHROM\\t%%POS\\t[%%END\\t%%GT\\t%%DR\\t%%DV\\t%%RR\\t%%RV]\\n' %(filename)s" % vars(self)):
			row = l.split()
			if row[3]!="1/1":continue
			results.append(row)
		return results
	def overlap_bed(self,bed_file):
		results = []
		bed = load_bed(bed_file,[1,2,3,4,5],4)
		calls = self.get_robust_calls()
		for call in calls:
			set_call_pos = set(range(int(call[1]),int(call[2])))
			for region in bed:
				if bed[region][0]!=call[0]: continue
				set_region_pos = set(range(int(bed[region][1]),int(bed[region][2])))
				intersect = set_call_pos.intersection(set_region_pos)
				if len(intersect)>1:
					results.append({"region":region,"start":min(intersect),"end":max(intersect)})
		print(results)
