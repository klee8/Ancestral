
## Merge MGLO and Progressive Cactus data

### DATA:
(MGLO files from Murray Cox, Cactus from David Winter)
+ cactus.maf : contains tree in first lines   
+ gene_order.dat : MGLO numbered genes in each tip genome
+ geneorder.out  : MGLO gene number order for each ancestral chr
+ gene_order.tree : MGLO tree
+ gff : gff files for each tip genome

### AIM:
+ read in cactus.maf blocks (each paragraph starting with 'a' is a multiple sequence alignment)
+ dump blocks less than 100bp
+ read in gene_order.dat
+ read in gff files for each tip genome
+ pair up gene order numbers with positional info and gene names for each species
+ identify genes in each cactus block and pair with MGLO gene number
+ see if blocks can be ordered by MGLO gene order in ancestral geneome (A6 in cactus, A0 in MGLO).

### Workflow:

#### Prep files for R

.maf files are a pain to read into R. Here we remove headers, empty lines and 
lines starting with 'a' (that demarcates the start of a new alignment block) and
set a new column at the start of each line that shows the block number (instead
of 's'). This will make it easy to read into R and use the 'group_by' tidyverse 
function to manipulate blocks individually
```
fmt_cactus_for_R.pl  data/cactus.maf data/cactus_fmtR.maf
```

MGLO files are also weird to get into R. Here we put each chromosome list of 
genes on one tab-delimited line with a column for genome name, a column for
chromosome number and a column with a comma-delimited list of gene-order numbers
```
fmt_MGLO_for_R.pl data/gene_order.dat data/gene_order_fmtR.dat
fmt_MGLO_for_R.pl data/geneorder.out data/geneorder_fmtR.out
```

#### merge info

```
merge_MGLO_cactus.R
```
