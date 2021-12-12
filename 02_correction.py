import os
import glob
from pathlib import Path
import snakemake
import gzip

contamination_threshold=2
workdir: config['workdir']
## 
sample_list = {}
with open(config['file_names_txt'],'r') as f:
    for line in f:
        items = line.strip().split("\t")
        sample_list[items[0]] = items[1:]

sample_P1 = {}
contaminated_list=config['checkm_list_info']
with open(contaminated_list,'r') as f:
    for line in f:
        items = line.strip().split("\t")
        contamintion=float(items[2])
        Circle=items[8]
        sample_name=items[0].split("_")[0]
        if contamintion > 2 and Circle == "N": 
            sample_P1[sample_name]="" ## Chimera candidate sample

SAMPLES = list(sample_list.keys())
SAMPLES_P1=list(sample_P1.keys())
SAMPLES_P2 = list(set(SAMPLES).difference(set(SAMPLES_P1)))

rule all:
    input:
        expand("{sample}/3.correction/{sample}.corrected.polished.fasta", sample = SAMPLES_P1),
        expand("{sample}/3.correction/{sample}.noerror.polished.fasta", sample = SAMPLES_P2)

rule get_contaminated_list:
    input:
        contaminated_list
    output:
        "{sample}/3.correction/bigtig_contaminated_seq"
    shell:
        """
        awk '$3>2 && $9 =="N"' {input} | grep {wildcards.sample}_ | cut -f1 > {output}
        """


rule bam_idx:
    input:
        '{some}.sort.bam'
    output:
        '{some}.sort.bam.bai'
    threads: 50
    shell:
        """
        samtools index -@ {threads} {input}
        """

rule bwa_align:
    #Align long reads to a fasta, storing the result in .bam format.
    input:
        '{sample}/2.polished/{sample}.polished.fa', 
        lambda wildcards: sample_list[wildcards.sample][1],
        lambda wildcards: sample_list[wildcards.sample][2]
    output:
        '{sample}/3.correction/bam_file/{sample}.short.sort.bam'
    threads: 50
    resources:
        time=2
    params:
        index_dir = '{sample}/3.correction/bwa_index/',
        index = '{sample}/3.correction/bwa_index/{sample}.bwa2_index',
        samtools_threads = 12
    shell:
        """
        if [ ! -d {params.index_dir} ]; then
            mkdir {params.index_dir}
        fi
        bwa-mem2 index -p {params.index} {input[0]} 
        bwa-mem2 mem -t {threads} {params.index} {input[1]} {input[2]} | samtools view -bS -@ {params.samtools_threads} - | samtools sort -@ 16 -m10g - -o {output}
        """

rule minimap:
    input:
        "{sample}/2.polished/{sample}.polished.fa",
        lambda wildcards: sample_list[wildcards.sample][0]
    output:
        '{sample}/3.correction/bam_file/{sample}.long.sort.bam'
    threads: 50
    shell:
        """
        minimap2 -t {threads} -aLQx map-ont --secondary=no --sam-hit-only {input} | samtools sort --threads {threads} -m3g - -o {output}
        """

rule seq_corr:
    input:
        rules.get_contaminated_list.output,
        "{sample}/3.correction/bam_file/{sample}.short.sort.bam",
        "{sample}/3.correction/bam_file/{sample}.short.sort.bam.bai",
        "{sample}/3.correction/bam_file/{sample}.long.sort.bam",
        "{sample}/3.correction/bam_file/{sample}.long.sort.bam.bai"
    output:
        "{sample}/3.correction/corrected_bigtigs/corrected.polished.fa"
    threads: 15
    params:
        outdir= "{sample}/3.correction/corrected_bigtigs"
    shell:
        """
        extract_contaminated_bigtigs_fasta.sh -s {input[0]} -sb {input[1]} -lb {input[3]} -o {params.outdir}
        """

rule final_part1:
    input:
        "{sample}/3.correction/corrected_bigtigs/corrected.polished.fa"
    output:
        "{sample}/3.correction/{sample}.corrected.polished.fasta"
    shell:
        """
        ln -s $PWD/{input} {output}
        """

rule final_part2:
    input:
        '{sample}/2.polished/{sample}.polished.fa'
    output:
        "{sample}/3.correction/{sample}.noerror.polished.fasta"
    shell:
        """
        ln -s $PWD/{input} {output}
        """ 
