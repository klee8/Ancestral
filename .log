  

AIM: Merge MLGO and Cactus Ancestral genomes

MLGO output of gene order 
Ancestral geneome A6 (root of tree) from cactus output 


NOTE: There seem to be different trees for cactus and mglo?
cactus:# hal (((Eam:1,Ebr:1)Anc3:1,Efe:1)Anc1:1,((Eel728:1,Eel732:1)Anc4:1,(Etp:1,(Ecl:1,Ety:1)Anc6:1)Anc5:1)Anc2:1)Anc0;
mglo:((Eam708,EfeFl1)A1,((Ebr1,(Eel728,Eel732)A2)A3,(Etp76,(Ety1756,Ecl1605)A4)A5)A6)A7;


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







