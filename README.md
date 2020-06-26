  

AIM: Reconstruct ancestral genomes

MLGO output of gene order 
Progressive Cactus from alignment data  

### Cactus data
#### export the HAL alignment as MAF referenced at the root
DATA from Amelie
hal2maf cactus.hal cactus.maf

### MGLO info
DATA from Murray Cox:
+ ortholog files
  1 The gene order in modern species is in ‘gene_order.txt’ 
  2	the phylogeny was fixed from standard phylogeny methods ‘gene_order.phy’
+ MGLO files
- geneorder.out
- gene_order.dat
geneorder.out file has lists of gene orders for ancestral genomes. 
ancestral genome A6 is the genome of interest

### Checked the distribution of E.clarkii chromosome 2 genes on Ancestral genome 6
 
```
perl check_mglo_output_overlaps.pl Ecl_gene_order.dat A6_geneorder.out > Ecl_chr_distribution_on_A6_genome.txt
```


```
gene_distribution_tip_v_anc.txt

```




### Checked the contiguous blocks of E.clarkii chromosome 2 on all tip and ancestral genomes
```
perl flag_gene_blocks.pl Ecl_chr2.txt Ecl_gene_order.dat MGLO_tagged.out
```







