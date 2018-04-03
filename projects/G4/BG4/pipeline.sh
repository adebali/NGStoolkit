cutadapt -f fastq -e 0.1 -q 20 -O 3 -a CTGTCTCTTATACACATCT $fq > /dev/stdout 2> $bname.cutadapt.txt \
| bwa mem -M genome.fa /dev/stdin \
| samtools view -S -u -F2304 -q 10 -L hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed - \
| samtools sort - ${bname}.hg19.tmp &&
java -Xmx5g -jar ~/bin/picard.jar MarkDuplicates I=${bname}.hg19.bam O=$bamclean/${bname}.hg19.clean.bam M=$logdir/$bname.md.txt