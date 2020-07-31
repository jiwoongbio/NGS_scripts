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
	my ($chromosome, $start, $end, $strand, $name) = split(/\t/, $line, -1);
	push(@{$nameChromosomeStartEndStrandListHash{$name}}, [$chromosome, $start, $end, $strand]);
}
close($reader);
foreach my $name (sort keys %nameChromosomeStartEndStrandListHash) {
	my @chromosomeStartEndStrandList = @{$nameChromosomeStartEndStrandListHash{$name}};
	next if(grep {$_->[1] eq '' || $_->[2] eq ''} @chromosomeStartEndStrandList);
	my $sequence = '';
	foreach(@chromosomeStartEndStrandList) {
		my ($chromosome, $start, $end, $strand) = @$_;
		if($strand eq '+') {
			$sequence = $sequence . uc($db->seq($chromosome, $start, $end));
		}
		if($strand eq '-') {
			$sequence = uc($db->seq($chromosome, $end, $start)) . $sequence;
		}
	}
	print ">$name\n";
	print "$sequence\n";
}
