#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'r' => \(my $reverse = ''),
	'c' => \(my $complementary = ''),
);
my $fastaFile = shift(@ARGV);
my $db = Bio::DB::Fasta->new($fastaFile);
if(defined(my $sequence = $db->seq(@ARGV))) {
	$sequence = uc($sequence);
	$sequence = reverse($sequence) if($reverse);
	$sequence =~ tr/ACGT/TGCA/ if($complementary);
	print "$sequence\n";
}
