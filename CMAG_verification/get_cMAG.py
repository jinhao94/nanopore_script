#!/usr/bin/env python3
import os,sys
from collections import defaultdict

aln   = sys.argv[1]
#########################################
# USAGE: python search_for_circ_contigs.py <CONTIGS.FASTA>
#
# This writes out a set of blasr commands
# specific to the supplied contigs that
# can be executed separately (might
# take awhile to run).
#
#
# Any blasr output files containing 
# alignments (i.e. non-zero file size)
# are possibly circularized and should
# be further investigated.
#########################################

# First preprocess the sequences so that we can align two portions of each sequence to each other
res_dict_start = defaultdict(list)
res_dict_end = defaultdict(list)
res_dict_tmp = defaultdict(list)

infile = open(aln, "r")

# R_start_pre = "NA"
# R_end_pre = "NA"

for line in infile:
	lines = line.split("\t")
	Q_ctg, Q_start, Q_end = lines[0:3]
	R_len, R_start, R_end = lines[6:9]
	#print(R_len, R_start, R_end)
	end_threshold = int(R_len) - 10000
	if int(R_start) < 10000:
		# res_dict[Q_ctg].append(line)
		if Q_ctg not in res_dict_start:
			res_dict_start[Q_ctg]=lines
			R_start_pre = int(R_start)
			res_dict_tmp[Q_ctg].append("start")
		else:
			if int(R_start) < R_start_pre:
				res_dict_start[Q_ctg]=lines
			else:
				continue
	elif int(R_end) > end_threshold:
		if Q_ctg not in res_dict_end:
			res_dict_end[Q_ctg] = lines
			R_end_pre = int(R_end)
			res_dict_tmp[Q_ctg].append("end")
		else:
			if int(R_end) > R_end_pre:
				res_dict_end[Q_ctg] = lines
			else:
				continue

hang_cig = "NA"
hang_len_prev = "NA"
for key in res_dict_tmp:
	if "start" in res_dict_tmp[key] and "end" in res_dict_tmp[key]:
		if int(res_dict_end[key][3]) - int(res_dict_start[key][3]) > 0 :
			Q_direction = "+"
		else:
			Q_direction = "-"
		if Q_direction == "+":
			gap_len = int(res_dict_end[key][2]) - int(res_dict_start[key][3])
		elif Q_direction == "-":
			gap_len = int(res_dict_end[key][3]) - int(res_dict_start[key][2])
		if hang_cig == "NA" :
			hang_cig = key
			hang_len_prev = int(gap_len)
		elif abs(hang_len_prev) > abs(int(gap_len)):
			hang_cig = key
			hang_len_prev = int(gap_len)


# print(res_dict_start[hang_cig], res_dict_end[hang_cig], gap_len)
print(aln, hang_len_prev, sep="\t")
#print(aln)



