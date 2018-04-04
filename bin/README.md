# NGS Toolkit Scripts

## FASTA modules

### Generate nucleotide content for each position from FASTA file

```
fa2kmerAbundanceTable.py -i input.fa -k 1 -o table.txt
```

`-k` indicates the length of the motif to count at each position. By default single nucleotide (A, T, C and G) is counted. It can be changes to any positive integer. For example `-k 2` can be used to list dinucleotides.

For using `ggplot2` in R, using 'melted' data is easier. To do that you can use the relevant script:

```
fa2kmerAbundanceMeltedData.py -i input.fa -o meltedTable.txt
```

for both scripts using `--percentage` will generate the fraction rather than the total raw counts.

