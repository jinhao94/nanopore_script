#!/usr/bin/sh

<<"!"
JinHao, 2020.01.13
Extract the fasta from bam file and reassembly the contigs
!

threads=8
function display_help() {
    echo "Extract the fastq from bam file which covered the  start and end position of largest conitgs "
    echo " "
    echo "Usage: extract_the_bigtigs_fasta.sh -b in.bam -s bigtigs.id -o bigtigs_terminal_reads.fq.gz "
    echo " "
    echo "   -lb, --long_bam       Long reads mapped sort.bam."
    echo "   -sb, --short_bam      Short reads mapped sort.bam."
    echo "   -s, --seq_list        Contaminated seq list."
    echo "   -o, --out_dir         Outdir."
    echo "   -t, --threads         Number of threads [default: 32]"
    echo "   -h, --help            Show this message."
    echo " "
    exit 1
}

if [ $# -eq 0 ];then
    echo "Please input parameters!"
    display_help
fi


function check_commond(){
    if [ $? -eq 0 ]; then
        echo "Environment loaded successfully"
    else
        echo "failed"
        exit
    fi
}

while [ "$1" != "" ]; do
    case $1 in
        -lb | --long_bam )           shift
                                long_bam=$1
                                ;;
        -sb | --short_bam )           shift
                                short_bam=$1
                                ;;
        -o | --out_dir )        shift
                                outdir=$1
                                ;;
        -s | --seq_list)        shift
                                seq_list=$1
                                ;;
        -t | --threads)         shift
                                threads=$1
                                ;;
        * )                     display_help
                                exit 1
    esac
    shift
done

## check the snakemake environment
#source activate snakemake
#check_commond

## check bam file
if [ ! -f $bam ]; then
    echo "Input bam is not existed."
    display_help
fi

# Get absolute path for snakemake
declare -a PATH_ARRAY # declare an indexed array variable
PATH_ARRAY[0]=$seq_list
PATH_ARRAY[1]=$long_bam
PATH_ARRAY[2]=$short_bam

n=0
for MY_PATH in "${PATH_ARRAY[@]}"; do
    if [[ $MY_PATH != /* ]]; then
        PATH_ARRAY[$n]=${PWD}/${MY_PATH}
    fi
    let "n = $n + 1"
done


## reassembly contaminated assemblies 
snakemake -s /gfsdata/gridengine/links/Seq_correction.py --directory $outdir --config contaminated_list=${PATH_ARRAY[0]} long_bam=${PATH_ARRAY[1]} short_bam=${PATH_ARRAY[2]} -p --cores $threads

#extract the short and long reads mapped on corrected assemblies
extract_r1=${outdir}/extracted_read/extract_short.r1.fq.gz
extract_r2=${outdir}/extracted_read/extract_short.r2.fq.gz
repaired_r1=${outdir}/extracted_read/repaired_short.r1.fq.gz
repaired_r2=${outdir}/extracted_read/repaired_short.r2.fq.gz
extract_long=${outdir}/extracted_read/extract_long.fq.gz



if [ ! -d ${outdir}/extracted_read ]; then
    mkdir -p ${outdir}/extracted_read
fi

## short reads extraction
seq_name=`cat $seq_list | xargs` 
samtools view -h -@ 6 $short_bam $seq_name | samtools fastq - -1 $extract_r1 -2 $extract_r2
/gfsdata/gridengine/software/bbmap/repair.sh in=$extract_r1 in2=$extract_r2 out=$repaired_r1 out2=$repaired_r2
## long reads extraction
samtools view -h -@ 6 $long_bam  $seq_name | samtools fastq - | pigz > $extract_long

## polish corrected assemblies
snakemake -s /gfsdata/gridengine/links/Polish.py --directory $outdir --config assembly=corrected.fasta sample_id=corrected assembly=corrected.fasta reads1=extracted_read/repaired_short.r1.fq.gz reads2=extracted_read/repaired_short.r1.fq.gz long_reads=extracted_read/extract_long.fq.gz --cores $threads -p
