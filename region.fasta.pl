#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;

my ($regionFile, $fastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($fastaFile);
my @nameList = ();
my %nameSequenceHash = ();
open(my $reader, $regionFile);
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome, $start, $end, $strand, $name) = split(/\t/, $line, 5);
	unless(defined($nameSequenceHash{$name})) {
		push(@nameList, $name);
		$nameSequenceHash{$name} = '';
	}
	$nameSequenceHash{$name} = $nameSequenceHash{$name} . uc($db->seq($chromosome, $start, $end)) if($strand eq '+');
	$nameSequenceHash{$name} = uc($db->seq($chromosome, $end, $start)) . $nameSequenceHash{$name} if($strand eq '-');
}
close($reader);
foreach my $name (@nameList) {
	print ">$name\n";
	print "$nameSequenceHash{$name}\n";
}
