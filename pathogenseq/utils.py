import re

indelre = re.compile("(\w)[\+\-](\d+)(\w+)")
def recode_indels(indels):
	#["C+5CGGGG","G-1C"]
	sorted_indels = sorted([x for x in indels],key=lambda y:len(y))
	largest_del = sorted([x for x in indels if "-" in x],key=lambda y:len(y))

	if len(largest_del)==0:
		leftseq = indelre.search(sorted_indels[-1]).group(1)
		rightseq = ""
	else:
		leftseq = indelre.search(largest_del[-1]).group(1)
		rightseq = indelre.search(largest_del[-1]).group(3)
	refseq = leftseq+rightseq
	recoded_indels = []
	for i in indels:
		if "-" in i:
			indel_len = int(indelre.search(i).group(2))
			iseq = leftseq+rightseq[:-(len(rightseq)-indel_len)]
		elif "+" in i:
			iseq = leftseq+indelre.search(i).group(3)+rightseq
		else:
			iseq = i+rightseq
		recoded_indels.append(iseq)
	return (refseq,recoded_indels)
