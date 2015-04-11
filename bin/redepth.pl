#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;


my %opts;
GetOptions(\%opts,"i:s","d:s","o:s","b:s","t:s","m:s","v:s","m:s","y:s");
my $expl = "If there's CNV in genome, you might use this script to remove depth bias effect from CNV. Each region has normalized copy number 1 for diploid chomosome status. For triploid, the copy number is 1.5 for example.";
my $usage=<<"USAGE";
        Program : $0
        Explanation : $expl
        Contact : deqiangs\@bcm.edu
        Usage : $0 [options]
                -i 	<str>	input BAM file name. 
                -o	<str>	output BAM file name.
		-b 	<str>	CNV region file name. It's a 4 column plain text file with colum 1-4 as chrom, start, end, copy number prediction.
                -m 	<str>	if your bam file RNAME has mate information '/1' or '/2', set it to 1. Otherwise 0. Default 1 for BSMAP alignment. 
                -v	<boleen>	verbose output	[0 or 1, default 0] 
                

USAGE

die $usage unless ($opts{"i"} and $opts{"b"} );

#=pod
my $inf = $opts{"i"} ;
my $outf = $opts{"o"};
my $cnvf = $opts{"b"};
my $m = (defined $opts{"m"}) ? $opts{"m"} : 1;







my %read = ();
open(IN, "bamToBed -i $inf | intersectBed -a $cnvf -b stdin -wa -wb | ") or die;
while(<IN>){
    chomp;
    my @f = split(/\t/, $_);
    my $name = $f[7];
    if($m == 1){
        $name = substr($name, 0, -2);
    }
    
    my $copyNum = $f[3];
    $read{$name} = $copyNum;
}
close(IN);

open(IN, "samtools view -h $inf | ") or die;
while(<IN>){
    chomp;
    my @f = split(/\t/, $_);
    my $name = $f[0];

    if(defined $read{$name}){
        my $copyNum = $read{$name};
        my $randNum = rand(1);
        if($randNum > 1/$copyNum){
            next; ##if $randNum ~~ (0, 1/$copyNum) -> keep;if $randNum ~~ (1/$copyNum, 1) -> delete; 
        }
    }
    print "$_\n";
}
close(IN);


