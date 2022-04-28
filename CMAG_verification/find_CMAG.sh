#!/bin/sh

#echo $#
function display_help() {
    echo -e "find CMAG"
    echo -e "Usage: find_CMAG.sh -g genome.fasta -b long_reads.sort.bam -o outfolder -t 8 \n"
    echo -e "outfolder: folder for outfile (Attention! If the outfolder already exits, the program will remove it.)"
    echo -e "threads: number of threads (default 8)"
    exit
}

if [ $# -lt 4 ]; then
    echo -e "Please input something\n"
    display_help
fi 

threads=8

while [ "$1" != "" ]; do
    case $1 in
        -g | --genome)             shift
                                   genome=$1
                                   ;;
        -b | --bam )               shift
                                   bam=$1
                                   ;;
        -o | --outfolder)          shift
                                   outfolder=$1
                                   ;;
        -t | --threads)            shift
                                   threads=$1
                                   ;;
        * )                        display_help
                                   exit 1
    esac
    shift
done


##
genome_prefix=`ls -d $genome | rush echo {%.}`
genome_reads=${outfolder}/${genome_prefix}.fq.gz
genome_reassembly=${outfolder}/${genome_prefix}.flye

contig_name=`seq_len $genome |  awk '{print $1}'`
genome_size=`seq_len $genome |  awk '{print $2}'`
let genome_end_cut=genome_size-50000



echo "samtools view -h $bam ${contig_name}:1-50000  ${contig_name}:${genome_end_cut}-${genome_size} | samtools fastq - | pigz > $genome_reads"
#flye --meta --plasmids --nano-raw $genome_reads --threads 16 --out-dir $genome_reassembly



