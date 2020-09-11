
# create bedfile for Eam from gff file 
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Eam708_Epichloe_amarillans_NFE702_38224169_v2.gff -o data/Eam.bed

# create bedfile for Anc0 (cactus root ancestor) from Eam bed file
halLiftover data/cactus.hal Eam data/Eam.bed Anc0 data/Anc0.bed

# extract Anc0 and Eam fasta files from hal file
hal2fasta data/cactus.hal Anc0 > data/Anc0.fasta
hal2fasta data/cactus.hal Eam > data/Eama702.fasta

# extract sequences from Anc0 fasta using bed co-ordinates
# reduce to max length transcript and then to chr:from-to
# for faidx, region file must be in the format chr:from-to
R --slave --vanilla < scripts/longest_transcript.R
samtools faidx data/Anc0.fasta --region-file data/Anc0_regions.txt --output data/Anc0.cds --mark-strand sign

