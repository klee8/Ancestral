# reformat MGLO file to parse into R
# change from fasta-style format to one line per chromosome. Each line should 
# be tab-delim, with first column as genome name, second column as chr name, 
# and third column with comma-separated genome order 

#!usr/bin/perl -w
use strict;

my $MGLO = $ARGV[0];
my $out = $ARGV[1];

open(MGLO, "<$MGLO") || die "ERROR, couldn't open $MGLO: $!";
open(OUT, ">$out") || "ERROR, couldn't open $out: $!";

my $genome = "";
my $chr = "";
my $seq = "";
while(<MGLO>) {
    chomp;
    if ($_ =~ /#/)  { next; }
    if ($_ =~ />/) { 
        $genome = $_; $genome =~ s/>//; 
        $chr = 1;
        $seq = $_;
        $seq =~ s/>$genome//;
    }
    elsif ($_ =~ /$/) { 
        $seq .= $_;
        $seq =~ s/\s+/,/g;
        $seq =~ s/\$//;
        print OUT "$genome\t$chr\t$seq\n"; 
        $seq = ""; $chr++;
    }
    else {
        $seq .= $_
    }
}

exit;
