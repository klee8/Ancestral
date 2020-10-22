# generate CHRONICle files from gff and protein.fna files for a genome


#!usr/bin/perl -w
use strict;

my $gff = $ARGV[0] || die "ERROR, need a gff file: $!";
my $name = $ARGV[1] || die "ERROR, requires short species name for output file: $!";
open(GFF, "<$gff") || die "ERROR, can't open $gff :$!";
open(OUT, ">$name.def") || die "ERROR, can't open $name.def: $!";

# add header line in after re-sorting the file later
#print OUT "type\tname\tchr\tstart\tend\tstrand\tsens\tIDg\/chr\tIDg\/all\tIDf/all\n";

# grab info for NAME.def file from gff
# Note that currently most genome annotations only have genes
# CHRONicle also uses tRNA, ncRNA, rRNA, repeatR, LTR, pseudoG, centromere
my $chr_gene_counter=1;
my $genome_gene_counter=1;
my $feature_counter=1;
my $current_chr=0;
my ($type, $name, $chr, $start, $end, $strand);
my ($sens, $IDg_chr, $IDg_all, $IDf_all);

while(<GFF>) {
    if ( $_=~ /\tgene\t/) {
        chomp;
        my @line = split("\t", $_);
        # skip tRNAs <<<<< these are missing in the fasta files
        if (abs($line[4]-$line[3]) < 150) { next; }   # tRNA usually under 100 but some genes here were longer
        # get gene name
        $line[8] =~ s/ID=//;
        $line[8] =~ s/;.*//;
        # set chromosome variable 
        my $chr_num = "";
        # check chr number format (some include underscores)
        if ($line[0]=~ /_/) { my @temp = split("_", $line[0]); $chr_num = $temp[1]; }
        if ($line[0]=~ /chr/) { $chr_num = $line[0]; $chr_num =~ s/chr//g; }
        # skip mitochondrial chrs (very few genes annotated on these anyway) 
        if (($line[0]=~ /mtDNA/) | ($line[0]=~ /_m/)  ) { next; }
        # reset chromosome counter on new chromosomes
        if ( $chr_num != $current_chr ) { 
            $current_chr = $chr_num;
            $chr_gene_counter = 1;
        }
        ($type, $name, $chr, $start, $end, $strand) = ($line[2], $line[8], sprintf("%03d", $chr_num), $line[3], $line[4], $line[6]);
        ($sens, $IDg_chr, $IDg_all, $IDf_all) = (0, sprintf("%05d", $chr_gene_counter), sprintf("%05d", $genome_gene_counter), sprintf("%05d", $genome_gene_counter));
        print OUT "$type\t$name\t$chr\t$start\t$end\t$strand\t$sens\t$IDg_chr\t$IDg_all\t$IDf_all\n";
        $chr_gene_counter++;
        $genome_gene_counter++; 
    }
#    note, have not been including other features (some, but not all of the genomes have a handful of tRNAs)
#    add block to include this info in the NAME.def file if you get it later   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#    elsif ( $_=~ /tRNA|ncRNA|rRNA|repeatR|LTR|pseudoG|centromere/ ) {
#       ($type, $name, $chr, $start, $end, $strand) = ($line[2], $line[8], sprintf("%03d", $chr_num), $line[3], $line[4], $line[6]);
#    	($sens, $IDg_chr, $IDg_all, $IDf_all) = (0, sprintf("%05d", $chr_gene_counter), sprintf("%05d", $genome_gene_counter), sprintf("%05d", $feature_counter));
#        print OUT "$type\t$name\t$chr\t$start\t$end\t$strand\t$sens\t$IDg_chr\t$IDg_all\t$IDf_all\n";
#       $feature_counter++;
#    }
}

close GFF; 
close OUT;



