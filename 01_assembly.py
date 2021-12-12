import os
import glob
from pathlib import Path
import snakemake
import gzip

sample_list = {}
with open(config['file_names_txt'],'r') as f:
    for line in f:
        items = line.strip().split("\t")
        sample_list[items[0]] = items[1:]



SAMPLES = list(sample_list.keys())

workdir: config['workdir']

rule all:
    input:
        expand("{sample}/2.polished/{sample}.polished.fa", sample = SAMPLES)
        

rule run_flye:
    input:
        lambda wildcards: sample_list[wildcards.sample][0]
    threads: 60
    params:
        outdir =  '{sample}/1.assembly'
    output:
        '{sample}/1.assembly/assembly.fasta'
    shell:
        """
            flye --meta -t {threads} --nano-raw {input} -g 300m -o {params.outdir}
        """


rule run_polish:
    input:
        rules.run_flye.output,
        lambda wildcards: sample_list[wildcards.sample]
    threads:50
    params:
        outdir = "{sample}/2.polished"
    output:
        "{sample}/2.polished/genome.nextpolish.fasta"
    shell:
        """
            nextpolish.sh -a {input[0]} -lr {input[1]} -r1 {input[2]} -r2 {input[3]} -o {params.outdir}
            nextPolish {wildcards.sample}/2.polished/run.cfg
        """

rule run_rename:
    input:
        rules.run_polish.output
    threads:30
    output:
        "{sample}/2.polished/{sample}.polished.fa"
    shell:
        """
        clean_nextpolish_seq_name.py {wildcards.sample}/2.polished/genome.nextpolish.fasta {output} {wildcards.sample}
        """

