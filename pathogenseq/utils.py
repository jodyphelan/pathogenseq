import re
from .fasta import *
indelre = re.compile("(\w)[\+\-](\d+)(\w+)")
gap_char = "-"

def recode_indels(indels):
	#["C+5CGGGG","G-1C"]
	sorted_indels = sorted([x for x in indels],key=lambda y:len(y))
	largest_del = sorted([x for x in indels if gap_char in x],key=lambda y:len(y))

	if len(largest_del)==0:
		leftseq = indelre.search(sorted_indels[-1]).group(1)
		rightseq = ""
	else:
		leftseq = indelre.search(largest_del[-1]).group(1)
		rightseq = indelre.search(largest_del[-1]).group(3)
	refseq = leftseq+rightseq
	recoded_indels = []
	for i in indels:
		if gap_char in i:
			indel_len = int(indelre.search(i).group(2))
			iseq = leftseq+rightseq[:-(len(rightseq)-indel_len)]
		elif "+" in i:
			iseq = leftseq+indelre.search(i).group(3)+rightseq
		else:
			iseq = i+rightseq
		recoded_indels.append(iseq)
	return (refseq,recoded_indels)

def variants2vcf(var_file,seq1_file,seq2_file,prefix,vcf_file):
	seq1_dict = fasta(seq1_file)
	seq2_dict = fasta(seq2_file)
	good_dp = 20
	realign_min_length = 5
	min_flank = 100

	seq1_chrom_i = 0
	seq1_pos_i = 1
	seq1_i = 2
	seq2_i = 3
	seq2_pos_i = 4
	seq2_chrom_i = 5
	gap_char = "-"

	del_lines = []
	ins_lines = []
	indel_line_set = set()
	tmp = []
	prev_type = None
	prev_pos = None
	prev_seq1_chrom = None
	lines = [l.rstrip().split() for l in open(var_file).readlines()]
	for i in range(len(lines)):
		#seq1_chromosome	4215484	G	C
		#Chroosome 212835	C	-

		row = lines[i]

		pos = int(row[seq1_pos_i])
		seq1_chrom = row[seq1_chrom_i]
		if row[seq1_i]==gap_char or row[seq2_i]==gap_char:
			indel_line_set.add(i)
			if prev_pos==None:
				tmp = [i]
			elif row[seq1_i]==gap_char:
				if pos==prev_pos:
					tmp.append(i)
				else:
					if prev_type=="ins":
						ins_lines.append(tmp)
					else:
						del_lines.append(tmp)
					tmp = [i]
			elif  row[seq2_i]==gap_char:
				if pos==prev_pos+1:
					tmp.append(i)
				else:
					if prev_type=="ins":
						ins_lines.append(tmp)
					else:
						del_lines.append(tmp)
					tmp = [i]
			prev_type = "ins" if row[seq1_i]==gap_char else "del"
			prev_pos = pos
			prev_seq1_chrom = seq1_chrom

	variants = defaultdict(dict)

	# for indel_pos in del_lines:
	# 	seq1_positions = [int(lines[i][seq1_pos_i]) for i in indel_pos]
	# 	seq2_positions = [int(lines[i][seq2_pos_i]) for i in indel_pos]
	# 	bases = [lines[i][seq1_i] for i in indel_pos]
	# 	seq1_chrom = lines[0][seq1_chrom_i]
	# 	query_chrom = lines[0][seq2_chrom_i]
	# 	indel_size = len(seq1_positions)
	# 	flank_size = indel_size if indel_size>min_flank else min_flank
	#
	# 	seq1_start,seq1_end = seq1_positions[0]-1,seq1_positions[-1]
	# 	seq2_start,seq2_end = seq2_positions[0],seq2_positions[-1]
	# 	seq1_left_flank,seq1_right_flank  =  seq1_start-flank_size,seq1_end+flank_size
	# 	seq2_left_flank,seq2_right_flank  =  seq2_start-flank_size,seq2_end+flank_size
	#
	# 	if indel_size>=realign_min_length:
	# 		print "-"*40
	# 		print "Anslysing deletion from %s to %s" % (seq1_start,seq1_end)
	# 		print "Extracting from %s %s:%s-%s" % (seq1_file,seq1_chrom,seq1_left_flank,seq1_right_flank)
	# 		print "Extracting from %s %s:%s-%s" % (seq2_file,query_chrom,seq2_left_flank,seq2_right_flank)
	# 		tmp_file_in = "%s.tmp.in.fa" % prefix
	# 		tmp_file_out = "%s.tmp.out.fa" % prefix
	# 		O = open(tmp_file_in,"w")
	# 		O.write(">seq1\n%s\n" % (seq1_dict.get_seq(seq1_chrom,seq1_left_flank,seq1_right_flank)))
	# 		O.write(">seq2\n%s\n" % (seq2_dict.get_seq(query_chrom,seq2_left_flank,seq2_right_flank)))
	# 		O.close()
	# 		muscle_align(tmp_file_in,tmp_file_out)
	# 		seq1_cnt = seq1_left_flank
	# 		seq2_cnt = seq2_left_flank
	# 		aln_dict = fasta(tmp_file_out)
	# 		tmp_seq1_bases = []
	# 		tmp_seq2_bases = []
	# 		tmp_seq1_positions = []
	# 		tmp_seq2_positions = []
	# 		for tmp_i in aln_dict.loop_pos("seq1"):
	# 			seq1_cnt = seq1_left_flank+tmp_i
	# 			seq2_cnt = seq2_left_flank+tmp_i
	# 			if aln_dict.fa_dict["seq1"][tmp_i] != aln_dict.fa_dict["seq2"][tmp_i]:
	# 				tmp_seq1_bases.append(aln_dict.fa_dict["seq1"][tmp_i])
	# 				tmp_seq2_bases.append(aln_dict.fa_dict["seq2"][tmp_i])
	# 				tmp_seq1_positions.append(seq1_cnt)
	# 				tmp_seq2_positions.append(seq2_cnt)
	# 		if tmp_seq1_bases!=bases:
	# 			print [lines[i][seq1_i] for i in indel_pos]
	# 			print tmp_seq1_bases
	# 			print tmp_seq2_bases
	# 			print tmp_seq1_positions
	# 			print tmp_seq2_positions
	# 			print first_seq1_pos-1
	# 			quit()
	#
	# 		bases = tmp_seq1_bases
	# 		seq1_positions = tmp_seq1_positions
	# 	first_seq1_pos = seq1_positions[0]
	#
	# 	switch = True
	# 	start_pos = first_seq1_pos-1
	# 	while switch:
	# 		n = seq1_dict.get_seq(seq1_chrom,start_pos)
	# 		end_base = bases[-1]
	# 		if n==end_base:
	# 			start_pos-=1
	# 			bases.insert(0,bases.pop())
	# 		else:
	# 			switch=False
	# 	alt_seq = seq1_dict.get_seq(seq1_chrom,start_pos)
	# 	ref_seq = alt_seq+"".join(bases)
	# 	variants[seq1_chrom][start_pos] = (ref_seq,alt_seq,"1/1",good_dp)
	# for indel_pos in ins_lines:
	# 	positions = [lines[i][seq1_pos_i] for i in indel_pos]
	# 	first_seq1_pos = int(positions[0])
	# 	seq1_chrom = lines[0][seq1_chrom_i]
	# 	indel_size = len(positions)
	# 	print "-"*40
	# 	print "Anslysing insertion from %s to %s" % (positions[0],positions[-1])
	# 	bases = [lines[i][seq2_i] for i in indel_pos]
	# 	switch = True
	# 	start_pos = first_seq1_pos
	# 	while switch:
	# 		n = seq1_dict.get_seq(seq1_chrom,start_pos)
	# 		end_base = bases[-1]
	# 		if n==end_base:
	# 			start_pos-=1
	# 			bases.insert(0,bases.pop())
	# 		else:
	# 			switch=False
	# 	ref_seq = seq1_dict.get_seq(seq1_chrom,start_pos)
	# 	alt_seq = ref_seq+"".join(bases)
	# 	variants[seq1_chrom][start_pos] = (ref_seq,alt_seq,"1/1",good_dp)

	for i in set(range(len(lines)))-indel_line_set:
		row = lines[i]
		pos,ref_seq,alt_seq = row[seq1_pos_i:seq1_pos_i+3]
		seq1_chrom = row[seq1_chrom_i]
		variants[seq1_chrom][int(pos)] = (ref_seq,alt_seq,"1/1",good_dp)


	OUT = open(vcf_file,"w")
	OUT.write("""##fileformat=VCFv4.1
##reference=/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa
##contig=<ID=seq1_chromosome,length=4411532>
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s
""" % prefix)
	npos = None
	nseq1_chrom = None
	for seq1_chrom in sorted(variants):
		for pos in sorted(variants[seq1_chrom]):
			if nseq1_chrom!=seq1_chrom:
				if pos!=1:
					npos=1
					OUT.write("%s\t%s\t.\t%s\t.\t.\t.\tEND=%s;MinDP=20\tGT:DP\t0/0:%s\n" % (seq1_chrom,npos,seq1_dict.get_seq(seq1_chrom,npos),pos-1,good_dp))
				else:
					pass
			else:
				if pos!=npos+1:
					OUT.write("%s\t%s\t.\t%s\t.\t.\t.\tEND=%s;MinDP=20\tGT:DP\t0/0:%s\n" % (nseq1_chrom,npos+1,seq1_dict.get_seq(seq1_chrom,npos),pos-1,good_dp))
			var = variants[seq1_chrom][pos]
			if pos==165608: log(var)
			if var[0]=="N" or var[0]=="n" or var[1]=="N" or var[1]=="n":
				if pos==165608: log(var)
				OUT.write("%s\t%s\t.\t%s\t%s\t255\t.\t.\tGT:DP\t%s:%s\n" % (seq1_chrom,pos,var[0],".","./.","."))
			else:
				OUT.write("%s\t%s\t.\t%s\t%s\t255\t.\t.\tGT:DP\t%s:%s\n" % (seq1_chrom,pos,var[0],var[1],var[2],var[3]))
			npos = pos
			nseq1_chrom = seq1_chrom


	OUT.close()




def muscle_align(filename,outfile):
	cmd = "muscle -in %s -out %s" % (filename,outfile)
	run_cmd(cmd)
