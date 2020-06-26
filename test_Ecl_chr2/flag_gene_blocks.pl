#!usr/bin/perl -w
use strict;

my $query_file = $ARGV[0] || "Ecl_chr2.txt";
my $target_file = $ARGV[1] || "Ecl_gene_order.dat";
my $outfile = $ARGV[2] || "MGLO_tagged.out";

# Assign Query chromosomes to a line in an array
open(QUERY, "<$query_file") || die "ERROR, couldn't open $query_file: $!";

my %query;
my $genome;
while(<QUERY>){
    chomp;
    if ($_=~ '>') {$genome = $_; $genome=~ s/>//; next;}
    $_=~ s/-//g;
    $_=~ s/\$//g;
    push (@{$query{$genome}}, $_);
}

# for each chromosome in target gene, flag any genes that are on the query chromsome
open(TARGET, "<$target_file") || die "ERROR, couldn't open $target_file: $!";
open(OUT, ">$outfile") || die "ERROR, couldn't open $outfile: $!";

my $target_chr=1;
print "Query\tTargetChr\tTargetChrGenes\tQueryChr\tQueryChrGenes\tOverlap\n";
while(<TARGET>){
    if ($_=~ ">") { print OUT $_; next;}
    chomp;
    $_=~ s/-//g;
    $_=~ s/\$//g;
    my $target = $_;
    my @target_genes = split(" ", $_);
    my $targetlength = @target_genes;
    # Query genome loop
    foreach my $genome (sort keys %query){
        my $chr;
        my $totalcount = 0;
        my $query_chr_num = 1;
        my $querylength;
        # Query chromosome loop
        foreach my $chrseq (@{$query{$genome}}) {
            my $chrcount = 0;
            my @query_genes=split(" ", $chrseq); 
            $querylength=@query_genes;
            # Query gene loop
            foreach my $ortho (@query_genes) {
                $ortho=~ s/ //g;
                if ( $ortho eq "") { next; }
                # check target
                if ( $target=~ /(^|\s)$ortho(\s|$)/ ) { 
                    $chrcount++; $totalcount++;
                    # flag target for output
                    $target=~ s/(^|\s)$ortho(\s|$)/ \*$ortho\* /g;
                }
            }        
            print "$genome\t$target_chr\t$targetlength\t$query_chr_num\t$querylength\t$chrcount\n";                
            # convert to binary (present/absent = 1/0)
            my @temp=split(" ", $target);
            my $binary= "";
            foreach my $gene (@temp) {
                if ($gene=~ /\*/) { $binary = $binary."1";}
                else { $binary = $binary."0";} 
            }
            print OUT "$binary\$\n";
            $query_chr_num++;
        }
        print "$genome\t$target_chr\t$targetlength\tall\t-\t$totalcount\n";
   }
   $target_chr++; 
}
