
### Identify overlaps of genes in a query set and target set of MLGO input or output.
```
perl check_mglo_output_overlaps.pl example/tip_file.txt example/anc_file.txt > overlaps.txt
```
Outputs:
QueryGenome	
QueryChr	
QueryChrLength	
TargetGenome	
TargetChr	
TargetChrLength	
Overlap

Note there is a line for each chromosome in the query file for each chromosome in the target file and then a summary line for all chromosomes in the query against that target

### Graph overlaps in R
```
R slave --vanilla < graph_distribution_query_on_target.R 
```
(Not set up for general use)

### Visualise where genes from a particular chromosome end up in a genome/genomes
This just outputs 0 for absence and 1 for presence of a gene from the query chromosome in the target file 
```
perl flag_gene_blocks.pl example/query_chr.txt example/anc_file.txt
```
