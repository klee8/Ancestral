# take protein fasta file and header .def file for CHRONicle
# create protein .prt file with headers including the same info as .def file

#!usr/bin/perl -w
use strict;


my $fastafile = $ARGV[0];
my $def_file = $ARGV[1];
my $species = $ARGV[2];

open(FASTA, "<$fastafile") || die "ERROR, couldn't open $fastafile: $!";
open(DEF, "<$def_file") || die "ERROR, couldn't open $def_file: $!";
open(OUT, ">$species.prt") || die "ERROR, couldn't open $species.prt: $!";



my %def;
while(<DEF>){
    # for genes only
    chomp;
    if ($_=~ /^gene/) {
        my @line = split("\t", $_);
        my $gene = $line[1];
        $gene =~ s/ //g;
#        $gene =~ s/Name.*//g;
        my $header = join("\t",@line[1..9]);
        $def{$gene}=$header;
    }
}

close DEF;


# loop through fasta
my $flag = 0;
while(<FASTA>) {
    my $gene;
    if ( ($_=~ /^>(.*)-T1/) || ($_=~ /^>(.*)\s/) )  {
        $gene = $1;
        $flag = 0; 
#        print "$gene: $def{$gene}\n";  
        if (exists($def{$gene})) {
            $_ =~ s/$gene-T1 //;
            $_ =~ s/$gene/$def{$gene}/;
            print OUT $_;
            $flag = 1;
            next;
        }    
    }
    if ($flag = 1) {
        print OUT $_;
    } 
}

close FASTA;
close OUT;
