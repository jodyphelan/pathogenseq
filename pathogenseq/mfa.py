from Bio import AlignIO
from .utils import *
from tqdm import tqdm

def mauve_call_variants(ref_file,query_file,prefix):
	aln_file = "%s.xmfa" % prefix
	cmd = "progressiveMauve --output=%s %s %s" % (aln_file,ref_file,query_file)
	run_cmd(cmd)
	ref_dict = fasta(ref_file)
	query_dict = fasta(query_file)
	variant_file = "%s.vars" % prefix
	O = open(variant_file,"w")
	alignment = AlignIO.parse(open(aln_file),"mauve")
	for seg in alignment:
		seq_names = []
		seqs = {}
		for record in seg:
			seqs[record.id] = record.seq
			seq_names.append(record.id)
		if len(seq_names)<2: continue
		ref_pos = 0
		alt_pos = 0
		for i in tqdm(range(seg.get_alignment_length())):
			bases = [seqs[s][i] for s in seq_names]
			ref_pos+=1
			alt_pos+=1
			if bases[0] != bases[1]:
				if bases[0]=="-":
					ref_pos-=1

				elif bases[1]=="-":
					alt_pos-=1

				chrom1 = ref_dict.get_chrom_by_pos(int(seq_names[0].split("/")[-1].split("-")[0])+2)
				chrom2 = query_dict.get_chrom_by_pos(int(seq_names[1].split("/")[-1].split("-")[0])+2)
				pos1 = int(seq_names[0].split("/")[-1].split("-")[0])+ref_pos
				pos2 = int(seq_names[1].split("/")[-1].split("-")[0])+alt_pos


				O.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1,pos1,bases[0],bases[1],pos2,chrom2))
	O.close()
	vcf_file = "%s.vcf" % prefix
	variants2vcf(variant_file,ref_file,query_file,prefix,vcf_file)
