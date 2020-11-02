### MCscan identifies synteny blocks and graphs them

Aim:
1 Compare ancestral genomes from MLGO and progressive cactus
2 Look at structural changes between E.clarkii -> E.festucae Fl1  -> E.typhina var poae
3 Look at how E clarkii chr 2 maps to each of the other genomes and how these chrs evolve across the tree

Notes:
Make sure that input gff and genome.fna files have the same chromosome names
anchors and last files affect downstream analysis. Clean out these files for every iteration
E clarkii chr 2 seems to map to chrs 2-7 in E clarkii (may need to filter out these genes)
can flip chromosomes by adding a minus in the seqid file or by switching the beginning and end positions in the layout (both end up with the name overlapping the end of the chr - fix in inkscape)
Chromosomes are not proportional across genomes
+ get length of each chr for each genome in the /data folder: for i in *.fna; do awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' < $i; done
+ tabulate in excel and try to get them in proportion to one another (there doesn't seem to be a simple way to do this. 'beginning' is at 0.1 and 'end' is at 0.8, some of this length may include the name also. Tried getting the proportion of the max expected length and then multiplying by the available space [guestimated at 0.6] and adding the 'start' 0.09])

To Do:
identify Ecl chr2 in other chromosomes - can these be filtered? what is actually matching up? add x to cds names in a 'Eclx.fna' genome for more direct comparison.


