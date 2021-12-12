#!/usr/bin/python
import os
import sys
from Bio import SeqIO


if len(sys.argv) != 4 :
    print("Usage: clean_nextpolish_seq_name.py genome.nextpolish.fasta clean_fasta sample_name")
    sys.exit()

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
sample_name = sys.argv[3]

if not os.path.exists(input_fasta) :
    print ("input fasta not exists.")
    sys.exit()

output_fasta_wt = open(output_fasta, "w")

for seq in SeqIO.parse(input_fasta, "fasta"):
    name = sys.argv[3] + "_" + "_".join(seq.id.split("_")[0:-1])
    out_seq = str(seq.seq).upper()
    output_fasta_wt.write(">" + name + "\n" + out_seq + "\n")
    # fasta_out.write_file(seq)
