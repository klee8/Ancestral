# check in MCscan which chromosomes line up with E clarkii chr2
# (previously used MLGO comparison data)

# RESULTS from below script with filtered synteny blocks:
# Etp - 1,2,3,6
# Ecl - NA (can't compare self)
# Ety - NA (not enough data)
# Eel728 - 2
# Eel732 - 1
# Ebr - 3
# Efe437 1,6
# Efl1 - 1,6
# Eam - 3,5

#!usr/bin/bash


# first run bed and cds setup for the following genomes
# Eam Ebr Efe437 Efl1 Eel728 Eel732 Ecl Ety Etp


########################   CHECK CHRS that map to Eclarkii chr 2
cd results_Ecl_chr2

rm *anchors*
rm *last*

# make a new copy of E clarkii files with changed name to compare it against itself
cp Ecl.bed Eclx.bed
cp Ecl.cds Eclx.cds
sed -i 's/mrna/mrnaX/g' Ecl.bed
sed -i 's/mrna/mrnaX/g' Ecl.cds

# calculate synteny blocks for comparison with Ecl chr2
python3 -m jcvi.compara.catalog ortholog Ecl Eam --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Ebr --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Efe437 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Efl1 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Eel728 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Eel732 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Eclx --cscore=.99 --no_strip_name     # renamed Ecl mrna to mrnaX in Eclx files in order to run it against itself 
python3 -m jcvi.compara.catalog ortholog Ecl Ety --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Etp --cscore=.99 --no_strip_name


# make new anchor files for karyotype
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Ecl.anchors Ecl.Ecl.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Eam.anchors Ecl.Eam.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Ebr.anchors Ecl.Ebr.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Efe437.anchors Ecl.Efe437.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Efl1.anchors Ecl.Efl1.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Eel728.anchors Ecl.Eel728.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Eel732.anchors Ecl.Eel732.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Ety.anchors Ecl.Ety.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Etp.anchors Ecl.Etp.anchors.new



## Eam
echo '''Ecl1605_2
Eam708_1,Eam708_2,Eam708_3,Eam708_4,Eam708_5,Eam708_6,Eam708_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .58,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .8,       0,      , Eam, top, Eam.bed
# edges
e, 0, 1, Ecl.Eam.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Eam.pdf


## Ebr
echo '''Ecl1605_2
Ebr1_1,Ebr1_2,Ebr1_3,Ebr1_4,Ebr1_5,Ebr1_6,Ebr1_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .78,       0,      , Ebr, top, Ebr.bed
# edges
e, 0, 1, Ecl.Ebr.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Ebr.pdf

## Ecl
echo '''Ecl1605_2
Ecl1605_1,Ecl1605_2,Ecl1605_3,Ecl1605_4,Ecl1605_5,Ecl1605_6,Ecl1605_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .22,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .8,       0,      , Eclx, top, Eclx.bed
# edges
e, 0, 1, Ecl.Eclx.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Eclx.pdf


## Efe437 -- Ecl chr 2 maps to Efe437 chr 1 and 6
echo '''Ecl1605_2
EfeE437_1,EfeE437_2,EfeE437_3,EfeE437_4,EfeE437_5,EfeE437_6,EfeE437_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .59,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .8,       0,      , Efe437, top, Efe437.bed
# edges
e, 0, 1, Ecl.Efe437.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Efe437.pdf


## Eel728
echo '''Ecl1605_2
Eel728_1,Eel728_2,Eel728_3,Eel728_4,Eel728_5,Eel728_6,Eel728_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .63,       0,      , Eel728, top, Eel728.bed
# edges
e, 0, 1, Ecl.Eel728.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Eel728.pdf

## Eel732
echo '''Ecl1605_2
Eel732_1,Eel732_2,Eel732_3,Eel732_4,Eel732_5,Eel732_6,Eel732_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .68,       0,      , Eel732, top, Eel732.bed
# edges
e, 0, 1, Ecl.Eel732.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Eel732.pdf

## Efl1
echo '''Ecl1605_2
EfeFl1_1,EfeFl1_2,EfeFl1_3,EfeFl1_4,EfeFl1_5,EfeFl1_6,EfeFl1_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .58,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .8,       0,      , Efl1, top, Efl1.bed
# edges
e, 0, 1, Ecl.Efl1.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Efl1.pdf

## Etp
echo '''Ecl1605_2
Etp76_1,Etp76_2,Etp76_3,Etp76_4,Etp76_5,Etp76_6,Etp76_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .29,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .8,       0,      , Etp, top, Etp.bed
# edges
e, 0, 1, Ecl.Etp.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Etp.pdf


## Ety
echo '''Ecl1605_2
Ety1756_1,Ety1756_2,Ety1756_3,Ety1756_4,Ety1756_5,Ety1756_6,Ety1756_7''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , Ecl, top, Ecl.bed
 .4,     .1,    .58,       0,      , Ety, top, Ety.bed
# edges
e, 0, 1, Ecl.Ety.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_proportional_chrs_Eclchr2.Ety.pdf

cd ..
