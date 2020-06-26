#!usr/bin/perl -w
use strict;

my $query_file = $ARGV[0] || "Ecl_gene_order.dat";
my $target_file = $ARGV[1] || "A6_geneorder.out";

# Assign Ecl chromosomes to a line in an array
open(QUERY, "<$query_file") || die "ERROR, couldn't open $query_file: $!";

my %query_genomes;
my $genome;
while(<QUERY>){
    chomp;
    if ($_=~ '>') {$genome = $_; $genome=~ s/>//; next;}
    $_=~ s/-//g;
    $_=~ s/$//g;
    push (@{$query_genomes{$genome}}, $_);
}

# for each chromosome in targetestral gene, count the number of hits of each gene to each Ecl chr
open(TARGET, "<$target_file") || die "ERROR, couldn't open $target_file: $!";
my $target_chr=1;
my $Tgenome="";
print "QueryGenome\tQueryChr\tQueryChrLength\tTargetGenome\tTargetChr\tTargetChrLength\tOverlap\n";
while(<TARGET>){
    chomp;
    if ($_=~ ">") {$Tgenome = $_; $Tgenome=~ s/>//; next;}
    $_=~ s/-//g;
    $_=~ s/$//g;
    my @target_genes = split(" ", $_);
    my $targetlength = @target_genes;
    foreach my $Qgenome (sort keys %query_genomes){
        my $chr;
        my $totalcount = 0;
        my $query_chr_num = 1;
        my $querylength;
        foreach my $chrseq (@{$query_genomes{$Qgenome}}) {
            my $chrcount = 0;
            my @temp=split(" ", $chrseq); 
            $querylength=@temp;
            foreach my $ortho (@target_genes) {
                $ortho=~ s/ //g;
                if ( $ortho eq "") { next; }
                if ( $chrseq=~ /(^|\s)$ortho(\s|$)/ ) { $chrcount++; $totalcount++;
                }
            }
            print "$Qgenome\t$query_chr_num\t$querylength\t$Tgenome\t$target_chr\t$targetlength\t$chrcount\n";                
            $query_chr_num++;
        }
        print "$Qgenome\tall\t-\t$Tgenome\t$target_chr\t$targetlength\t$totalcount\n";
   }
   $target_chr++; 
}
