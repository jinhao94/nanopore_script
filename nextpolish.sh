#!/usr/bin/sh

#echo $#
function display_help() {
    echo -e "jinh 20200421"
    echo -e "Prior to utilizing this tool, it's essential to first install the NextPolish software. Following installation, please adjust the 'script_path' variable in the settings to reflect the absolute path where NextPolish is located."
    echo -e "The read path must be absolute path"
    echo -e "Usage: nextpolish.sh -r1 sceondary_assembly.fasta -r2 primary_assembly.fasta -lr -o  outfolder \n"
    echo -e "-r1 : short r1 reads "
    echo -e "-r2 : short r2 reads "
    echo -e "-lr : long reads "
    echo -e "-a : assembly"
    echo -e "-o : folder for outfile (Attention! If the outfolder already exits, the program will remove it.)"
    # echo -e "threads: number of threads (default 8)"
    exit
}


if [ "$#" == 0 ]; then
    display_help
fi

while [ "$1" != "" ]; do
    case $1 in
        -r1 | --short_r1)          shift
                                   short_r1=$1
                                   ;;
        -r2 | --short_r2)          shift
                                   short_r2=$1
                                   ;;
        -lr | --longreads)         shift
                                   longreads=$1
                                   ;;
        -a | --assembly)           shift
                                   assembly=$1
                                   ;;
        -o | --outfolder)          shift
                                   outfolder=$1
                                   ;;
        * )                        display_help
                                   exit 1
    esac
    shift
done

script_path='/nvmessdnode3/opt/software/NextPolish'
wkdir=./

if [ ! -d $outfolder ];then
        mkdir $outfolder
else
        echo "Out folder already exist, so it will be removed!"
        rm -rf $outfolder
        mkdir $outfolder
fi

# mkdir ${outfolder}/01_short_reads_rawdata ${outfolder}/02_long_reads_rawdata

function path_convert() {
    echo $1 | perl -e 'while(<STDIN>){chomp; if($_=~/^\//){print "$_"}else{print "@ARGV[0]/$_"} } ' $PWD
    # echo $1
}


r1p=`path_convert $short_r1`
r2p=`path_convert $short_r2`
lp=`path_convert $longreads`
genome=`path_convert $assembly`

sgc=./sgc.config
lgc=./lgc.config

echo -e "$r1p\n$r2p" > ${outfolder}/$sgc
echo $lp > ${outfolder}/$lgc
# echo $r1p $r2p $lp
perl -e 'open I, @ARGV[0]; while(<I>){chomp; if($_=~/^genome =/){print "genome = @ARGV[1]\n"; next}; if($_=~/^workdir/){print "workdir = @ARGV[2]\n"; next}; if($_=~/^sgs_fofn/){print "sgs_fofn = @ARGV[3]\n"; next}; if($_=~/lgs_fofn =/){print "lgs_fofn = @ARGV[4]\n"; next}; print "$_\n" }' ${script_path}/run.cfg $genome $wkdir $sgc $lgc > ${outfolder}/run.cfg
