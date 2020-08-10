# number blocks in maf file to parse into R
#!usr/bin/perl -w
use strict;

my $cactus = $ARGV[0];

open(CACTUS, "<$cactus") || die "ERROR, couldn't open $cactus: $!";
open(OUT, ">cactus.fmtR.maf") || "ERROR, couldn't open outfile: $!";

my $counter = 0;
while(<CACTUS>) {
    if ($_ =~ /^a/) { $counter ++; next;}
    if ($_ =~ /^\s/) { next; }
    if ($_ =~ /#/)  { next; }
    print OUT "$counter\t$_";
}

exit;
