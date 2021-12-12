import os
import glob
from pathlib import Path
import snakemake
import gzip

#workdir: config['workdir']
seq_list = {}
with open(config['contaminated_list'],'r') as f:
    for line in f:
        items = line.strip().split("\t")
        seq_list[items[0]] = items[1:]

SEQ = list(seq_list.keys())
print(SEQ)
long_bam=config['long_bam']
short_bam=config['short_bam']

rule all:
    input:
        expand("corrected.fasta")


rule extract_long_fasta:
    input:
        long_bam
    output:
        "seq_reads/{seq}.fasta.gz"
    shell:
        """
        samtools view -h {input} {wildcards.seq} | samtools fasta | pigz > {output}
        """

rule re_assembly:
    input:
        rules.extract_long_fasta.output
    output:
        "seq_reassembly/{seq}/assembly.fasta"
    threads: 12
    params:
        outdir='seq_reassembly/{seq}'
    shell:
        """
        flye --meta -m 3000 -t {threads} --nano-raw {input} -o {params.outdir}
        """


rule rename_assembly:
    input:
        rules.re_assembly.output
    output:
        'seq_reassembly.links/{seq}.fa'
    run:
        import sys, os
        from Bio import SeqIO
        tmp = 0
        outfile = open(output[0], 'w')
        for seq_record in SeqIO.parse(input[0], "fasta"):
            seq_name=wildcards.seq
            tmp += 1
            outfile.write(">" + seq_name + "_corr_" + str(tmp) + "\n" + str(seq_record.seq) + "\n")

rule final:
    input:
        expand("seq_reassembly.links/{seq}.fa", seq = SEQ)
    output:
        "corrected.fasta"
    shell:
        """
        cat seq_reassembly.links/*fa > corrected.fasta
        """


