import os
import glob
from pathlib import Path
import snakemake
import gzip

#workdir: config['workdir']
sample_id = config['sample_id']
assembly = config['assembly']
reads1 = config['reads1']
reads2 = config['reads2']
long_reads = config['long_reads']

long_read_polish_round=4

rule all:
    input:
        expand("{sample}.polished.fa", sample = sample_id)


def get_long_polish_input(wildcards):
    #this method will choose a  iteration for the unpolished assembly, depending on the value of the iteration wildcard
    if int(wildcards.iterations) == 1:
        sample_name = assembly
        return(sample_name, sample_name + ".long.sort.bam" , sample_name + ".long.sort.bam.bai", sample_name + ".fai")
    else:
        result = "long_polish/long_{iterations}/{sample}_long_{iterations}.fa".format(sample = wildcards.sample, iterations = str(int(wildcards.iterations) - 1))
        return(result, result + ".long.sort.bam", result + ".long.sort.bam.bai", result + ".fai")

rule long_read_polish:
    input:
        get_long_polish_input
    output:
        'long_polish/long_{iterations}/{sample}_long_{iterations}.fa'
    shell:
        """
        bam_list=long_polish/long_{wildcards.iterations}/input.long.sort.bam.list
        echo {input[1]} > $bam_list
        python /gfsdata/gridengine/software/NextPolish/lib/nextpolish2.py -r ont -g {input[0]} -p {threads} -l $bam_list > {output}
        """


def get_short_t1_input(wildcards):
    assembly= "long_polish/long_{iterations}/{sample}_long_{iterations}.fa".format(sample=sample_id, iterations=long_read_polish_round)
    return(assembly, assembly + ".short.sort.bam", assembly + ".short.sort.bam.bai", assembly + ".fai")

rule run_short_polish_task_1:
    input:
        get_short_t1_input
    output:
        'short_polish/{sample}_short_task1.fa'
    threads: 16
    shell:
        """
        python /gfsdata/gridengine/software/NextPolish/lib/nextpolish1.py -ploidy 1 -g {input[0]} -t 1 -p {threads} -s {input[1]} > {output} 
        """


rule run_short_polish_task_2:
    input:
        'short_polish/{sample}_short_task1.fa',
        'short_polish/{sample}_short_task1.fa.short.sort.bam',
        'short_polish/{sample}_short_task1.fa.short.sort.bam.bai'
    output:
        'short_polish/{sample}_short_task2.fa'
    threads: 16
    shell:
        """
        #mkdir short_polish_task2
        /usr/bin/python /gfsdata/gridengine/software/NextPolish/lib/nextpolish1.py -ploidy 1 -g {input[0]} -t 2 -p 1 -s {input[1]} > {output}
        """

rule final:
    input:
        rules.run_short_polish_task_2.output
    output:
        "{sample}.polished.fa"
    shell:
        """
        ln -s $PWD/{input} {output}
        """


rule align_long_reads:
    input:
        '{some}',
        long_reads
    output:
        '{some}.long.sort.bam',
        '{some}.long.sort.bam.bai'
    threads: 30
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input} | samtools sort - -m 2g --threads {threads} -o {output[0]};
        samtools index -@ 6 {output[0]}
        """


rule align_short_reads:
    #Align short reads to the unpolished or longread-polished assembly, depending on the output of
    input:
        '{some}',
        reads1,
        reads2
    output: 
        "{some}.short.sort.bam", 
        "{some}.short.sort.bam.bai"
    threads:30
    shell:
        """
        bwa-mem2 index {input[0]}; 
        bwa-mem2 mem -t {threads} {input[0]} {input[1]} {input[2]} | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3 - - | samtools sort --threads {threads} - | samtools markdup --threads 5 -r - {output[0]}
        samtools index {output[0]}
        """


rule samtools_idx:
    input:
        '{some}'
    output:
        '{some}.fai'
    threads: 12
    shell:
        """
        samtools faidx {input} 
        """
