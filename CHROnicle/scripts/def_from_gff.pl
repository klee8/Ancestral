# generate CHRONICle files from gff and protein.fna files for a genome


#!usr/bin/perl -w
use strict;

my $gff = $ARGV[0] || die "ERROR, need a gff file: $!";
my $name = $ARGV[1] || die "ERROR, requires short species name for output file: $!";
open(GFF, "<$gff") || die "ERROR, can't open $gff :$!";
open(OUT, ">$name.def") || die "ERROR, can't open $name.def: $!";

# grab info for NAME.def file from gff
# Note that currently the annotation only has genes
# CHRONicle also uses tRNA, ncRNA, rRNA, repeatR, LTR, pseudoG, centromere
my $chr_counter=1;
my $genome_counter=1;
my $feature_counter=1;
my $current_chr=1;
my ($type, $name, $chr, $start, $end, $strand);
my ($sens, $IDg_chr, $IDg_all, $IDf_all);
while(<GFF>) {
    if ( $_=~ /\tgene\t/) {
        chomp;
        my @line = split("\t", $_);
        $line[8] =~ s/ID=//;
        $line[8] =~ s/;.*//; 
        my $chr_num = "";
        if ($line[0]=~ /_/) { my @temp = split("_", $line[0]); $chr_num = $temp[1]; }
        if ($line[0]=~ /chr/) { $chr_num = $line[0]; $chr_num =~ s/chr//g; }
        if ($line[0]=~ /mtDNA/) { $chr_num = 8; }
        # reset chromosome counter on new chromosomes
        if ( $chr_num != $current_chr ) { 
            $current_chr = $chr_num;
            $chr_counter = 1;
        }
        ($type, $name, $chr, $start, $end, $strand) = ($line[2], $line[8], $chr_num, $line[3], $line[4], $line[6]);
        ($sens, $IDg_chr, $IDg_all, $IDf_all) = (0, $chr_counter, $genome_counter, $feature_counter);
        $chr_counter++;
        $genome_counter++; 
        $feature_counter++;
    }
    elsif ( $_=~ /tRNA|ncRNA|rRNA|repeatR|LTR|pseudoG|centromere/ ) {

    # add block to include this info in the NAME.def file if you get it later   <<<<<<<<<<<<<<<<
        $feature_counter++;
    }
    else { next; }
    print OUT "$type\t$name\t$chr\t$start\t$end\t$strand\t$sens\t$IDg_chr\t$IDg_all\t$IDf_all\n";
}

close GFF; 
close OUT;



