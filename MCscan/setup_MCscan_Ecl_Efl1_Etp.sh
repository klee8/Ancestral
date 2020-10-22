#!usr/bin/bash

## rename chrs to match between gff and genome.fna files
#sed -i 's/Ecl_1605_22/Ecl1605/g' data/Epichloe_clarkii_Hl.gff3
#sed -i 's/chr/EfeFl1/g' data/EfFl1_v3.1.gff


#######################################################
#            MAKE CACTUS ROOT BED AND CDS FILES       #
#######################################################
echo "creating bed and cds files.............."

# create bedfiles from gff files
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_clarkii_Hl.gff3 -o data/Ecl.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/EfFl1_v3.1.gff -o data/Efl1.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_typhina_poae_NFe76.gff -o data/Etp.bed

# extract sequences from genome fasta using bed co-ordinates
# reduce to max length transcript and then to chr:from-to
# for faidx, region file must be in the format chr:from-to
for i in Ecl Efl1 Etp;
do
    Rscript --vanilla  scripts/longest_transcript.R data/$i.bed
done

# replace old bedfile with longest transcript bedfile
for i in  Ecl Efl1 Etp
do
    mv data/$i.maxlen.bed data/$i.bed
done

# get cds data for each file using samtools
samtools faidx data/Ecl1605_22_Epichloe_clarkii_1605_22_45646793_v1.fna --region-file data/Ecl.regions.txt --output data/Ecl.cds  --mark-strand sign
samtools faidx data/EfeFl1_Epichloe_festucae_Fl1_35023690_v2.fna --region-file data/Efl1.regions.txt --output data/Efl1.cds --mark-strand sign
samtools faidx data/Etp76_Epichloe_typhina_var_poae_NFe76_38327242_v1.fna --region-file data/Etp.regions.txt --output data/Etp.cds --mark-strand sign

# note [faidx] Truncated sequence: Ecl1605_3:7592542-7597818 for E clarkii

# pop mRNA_# data back into cds headers from lookup file
for i in Ecl Efl1 Etp; 
do
    # linearise cds file
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < data/$i.cds > data/$i.lin.cds
    # replace with region lookup file info
    sed 's/mrna/>mrna/g' data/$i.regions_lookup.txt > data/$i.regions_headers.txt
    cut -f 2 data/$i.lin.cds > data/$i.seqs.cds
    paste data/$i.regions_headers.txt data/$i.seqs.cds > data/$i.cds
    # use perl to change the second instance of tab to a carrige return in each line
    perl -pe 's{\t}{++$n % 2 ? $& : "\n"}ge' data/$i.cds > temp
    mv temp data/$i.cds
    # re-format to 60 chars per line
    cp data/$i.cds .
    R --slave --vanilla < scripts/format_60char_fasta.R --args $i.cds
    mv fmt_$i.cds data/$i.cds
    rm $i.cds
    # MCscan takes the first number in a chromosome name as the chromosome number. 
    # change chromosome names as appropriate to accommodate this
done


########################################################
#              USE MCscan
#########################################################

echo "create and populate results directory ................."
mkdir results
cd results
for i in Ecl Efl1 Etp
do
    cp ../data/$i.bed .
    cp ../data/$i.cds .
done

# format CDS files (to transcript name only)
#for i in Ecl Efl1 Etp
#do
#    python -m jcvi.formats.fasta format $i.cds > tmp.cds
#    mv tmp.cds $i.cds
#done

# synteny blocks
echo "calculate synteny blocks ................."
python -m jcvi.compara.catalog ortholog Ecl Efl1 --no_strip_names
python -m jcvi.compara.catalog ortholog Efl1 Etp --no_strip_names

# filter synteny blocks
#echo "calculate filtered synteny blocks............."
rm Ecl.Efl1.last.filtered
rm Efl1.Etp.last.filtered
python3 -m jcvi.compara.catalog ortholog Ecl Efl1 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Efl1 Etp --cscore=.99 --no_strip_name

# compare synteny reciprocal depth of coverage <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Doesn't work for this data set?????
echo "compare synteny reciprocal depth of coverage.............."
python -m jcvi.compara.synteny depth --histogram Ecl.Efl1.anchors
python -m jcvi.compara.synteny depth --histogram Efl1.Etp.anchors

# make dotplot (already done by filter synteny block command)
echo "make dotplot ............"
python -m jcvi.graphics.dotplot Ecl.Efl1.anchors
python -m jcvi.graphics.dotplot Efl1.Etp.anchors


##### manually make seqids and layout plots  ############

### seqids < list of comma separated chr names, one spp per line
#chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19
#Pp01,Pp02,Pp03,Pp04,Pp05,Pp06,Pp07,Pp08

echo '''Ecl1605_1,Ecl1605_2,Ecl1605_3,Ecl1605_4,Ecl1605_5,Ecl1605_6,Ecl1605_7
EfeFl1_1,EfeFl1_2,EfeFl1_3,EfeFl1_4,EfeFl1_5,EfeFl1_6,EfeFl1_7
Etp76_1,Etp76_2,Etp76_3,Etp76_4,Etp76_5,Etp76_6,Etp76_7'''> seqids

### layout < example graph layout
## y, xstart, xend, rotation, color, label, va,  bed
# .6,     .1,    .8,       0,      , Grape, top, grape.bed
# .4,     .1,    .8,       0,      , Peach, top, peach.bed
## edges
#e, 0, 1, grape.peach.anchors.simple

echo '''# y, xstart, xend, rotation, color, label, va,  bed
.7,     .1,    .8,       0,       , Ecl, top, Ecl.bed
.5,     .1,    .8,       0,       , Efl1, top, Efl1.bed
.3,     .1,    .8,       0,       , Etp, top, Etp.bed
# edges
e, 0, 1, Ecl.Efl1.anchors.simple
e, 1, 2, Efl1.Etp.anchors.simple''' > layout


# make karyotype plot
#echo "make karyotype plot..............."
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Efl1.anchors Ecl.Efl1.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Efl1.Etp.anchors Efl1.Etp.anchors.new
python3 -m jcvi.graphics.karyotype seqids layout

cd ..


####### SORTED figures by hand
#Note: removed Chr18,Chr87,Chr185,Chr121,Chr130,Chr197,Chr175,Chr213 from ssid list for MLGO chr 5 as they were duplicated from MLGO chr 1
 
