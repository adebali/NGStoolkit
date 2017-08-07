# Single End vs Paired End Pipeline Differences

## Samtools view
PE: '-bf 0x2'

SE: '-bF 0x4'

## bam2bed _bedtools bamtobed_
PE: '-bedpe -mate1'

SE: No flags
