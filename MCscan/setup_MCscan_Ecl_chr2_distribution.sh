#!usr/bin/bash

#mkdir results_chr2
cd results_chr2
rm *anchors*
rm *last*

for i in Eam Ebr Efe437 Efl1 Eel728 Eel732 Ecl Ety Etp
do
    cp ../data/$i.bed .
    cp ../data/$i.cds .
done

 format CDS files (to transcript name only)
for i in Eam Ebr Efl1 Eel728 Eel732 Ecl Ety Etp
do
    python -m jcvi.formats.fasta format $i.cds tmp.cds
    mv tmp.cds $i.cds
done

# calculated filtered synteny blocks
# remove intermediate files from prior iterations
rm *anchors*
rm *last*
#echo "calculate filtered synteny blocks............."
python3 -m jcvi.compara.catalog ortholog Eam Ebr --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ebr Efe437 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Efe437 Efl1 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Efl1 Eel728 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Eel728 Eel732 --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Eel732 Ecl --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Ety --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ety Etp --cscore=.99 --no_strip_name
python3 -m jcvi.compara.catalog ortholog Ecl Efe437 --cscore=.99 --no_strip_name

# compare synteny reciprocal depth of coverage 
#echo "compare synteny reciprocal depth of coverage.............."
python -m jcvi.compara.synteny depth --histogram Eam.Ebr.anchors
python -m jcvi.compara.synteny depth --histogram Ebr.Efl1.anchors
python -m jcvi.compara.synteny depth --histogram Efl1.Eel728.anchors
python -m jcvi.compara.synteny depth --histogram Eel728.Eel732.anchors
python -m jcvi.compara.synteny depth --histogram Eel732.Ecl.anchors
python -m jcvi.compara.synteny depth --histogram Ecl.Ety.anchors
python -m jcvi.compara.synteny depth --histogram Ety.Etp.anchors

# make dotplot (already done by filter synteny block command)
#echo "make dotplot ............"
#python -m jcvi.graphics.dotplot Eam.Ebr.anchors
#python -m jcvi.graphics.dotplot Ebr.Efl1.anchors
#python -m jcvi.graphics.dotplot Efl1.Eel728.anchors
#python -m jcvi.graphics.dotplot Eel728.Eel732.anchors
#python -m jcvi.graphics.dotplot Eel732.Ecl.anchors
#python -m jcvi.graphics.dotplot Ecl.Ety.anchors
#python -m jcvi.graphics.dotplot Ety.Etp.anchors

# make new anchor files for karyotype
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Eam.Ebr.anchors Eam.Ebr.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ebr.Efe437.anchors Ebr.Efe437.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Efe437.Efl1.anchors Efe437.Efl1.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Efl1.Eel728.anchors Efl1.Eel728.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Eel728.Eel732.anchors Eel728.Eel732.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Eel732.Ecl.anchors Eel732.Ecl.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Ety.anchors Ecl.Ety.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ety.Etp.anchors Ety.Etp.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Efe437.anchors Ecl.Efe437.anchors.new


######## manually make seqids and layout plots  ############

# MCscan can only map three genomes at a time
# divided data into 4 overlapping groups of three in the same order as the tree

# seqid
# < list of comma separated chr names, one spp per line
# < list of chrs that map to E. clarkii chr 2 (from Ecl_chr2 mapping)

# layout
# < Changed xend to be in proportion to sum of Etyp chr lengths used below 
# < see end_positions_for_proportional_MCscan_layout.xlxs file for details

# files for Eam Ebr Efe437
echo '''Eam708_3,Eam708_5
Ebr1_3
EfeE437_1,EfeE437_2''' > seqids

# layout file. Changed xend to be in proportion to sum of Etyp chr lengths used below 
echo '''# y, xstart, xend, rotation, color, label, va,  bed
.7,     .1,    .32,       0,       , Eam, top, Eam.bed
.5,     .1,    .25,       0,       , Ebr, top, Ebr.bed
.3,     .1,    .32,       0,       , Efe437, top, Efe437.bed
# edges
e, 0, 1, Eam.Ebr.anchors.simple
e, 1, 2, Ebr.Efe437.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_Eam.Ebr.Efe437.pdf


# files for Efe437 Efl1 Eel728
echo '''EfeE437_1,EfeE437_2
EfeFl1_1,EfeFl1_2,EfeFl1_6
Eel728_2''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
.7,     .1,    .32,       0,       , Efe437, top, Efe437.bed
.5,     .1,    .32,       0,       , Efl1, top, Efl1.bed
.3,     .1,    .21,       0,       , Eel728, top, Eel728.bed
# edges
e, 0, 1, Efe437.Efl1.anchors.simple
e, 1, 2, Efl1.Eel728.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_Efe437.Efl1.Eel728.pdf


# files for Eel728 Eel732 Ecl
echo '''Eel728_2
Eel732_1
Ecl1605_2''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
.7,     .1,    .21,       0,       , Eel728, top, Eel728.bed
.5,     .1,    .22,       0,       , Eel732, top, Eel732.bed
.3,     .1,    .26,       0,       , Ecl, top, Ecl.bed
# edges
e, 0, 1, Eel728.Eel732.anchors.simple
e, 1, 2, Eel732.Ecl.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_Eel728.Eel732.Ecl.pdf


# files for Ecl Ety Etp
echo '''Ecl1605_2
Ety1756_2
Etp76_1,Etp76_2,Etp76_3,Etp76_4,Etp76_5,Etp76_6''' > seqids

echo '''# y, xstart, xend, rotation, color, label, va,  bed
.7,     .1,    .25,       0,       , Ecl, top, Ecl.bed
.5,     .1,    .2,       0,       , Ety, top, Ety.bed
.3,     .1,    .69,       0,       , Etp, top, Etp.bed
# edges
e, 0, 1, Ecl.Ety.anchors.simple
e, 1, 2, Ety.Etp.anchors.simple''' > layout

python3 -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_Ecl.Ety.Etp.pdf


cd ..


