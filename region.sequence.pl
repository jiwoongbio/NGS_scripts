#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;

my ($regionFile, $fastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($fastaFile);
my %nameChromosomeStartEndStrandListHash = ();
open(my $reader, $regionFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my ($chromosome, $start, $end, $strand, $name) = @tokenList;
	my $sequence = '';
	$sequence = uc($db->seq($chromosome, $start, $end)) if($strand eq '+');
	$sequence = uc($db->seq($chromosome, $end, $start)) if($strand eq '-');
	print join("\t", @tokenList, $sequence), "\n";
}
close($reader);
