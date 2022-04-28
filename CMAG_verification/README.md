ls -d CMAG.bam/*bam | parallel --plus -j 5 sh find_CMAG.sh -g CMAG.fas/{/..} -b {} -o CMAG.fastq > run_extract.sh
samtools view -h A102.CMAG_1.fa.sort.bam A102_contig_1238:1-50000  A102_contig_1238:2336017-2386017 | samtools fastq - | pigz > A102.CMAG_1.gap.reads.fq.gz
flye  --meta --subassemblies A102.CMAG_1.fa -o A102.CMAG_1.gap.ass.meta/ -t 32
ls -d CMAG.flye/*/assembly.fasta | rush -d / --dry-run 'minimap2 -PD -k19 -w19 -m200 -t8 CMAG.fas/{2}.fa {} > CMAG.aln/{2}.aln ' | sed 's/\.flye.fa/.fa/g'
ls -d *aln | parallel get_cMAG.py {}
