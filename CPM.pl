#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'H' => \(my $header = ''),
);
my ($geneCountFile) = @ARGV;
die "Gene count file is not defined.\n" unless(defined($geneCountFile));
die "$geneCountFile is not readable.\n" unless(-r $geneCountFile);
my @geneCountListList = ();
my @totalCountList = ();
{
	open(my $reader, $geneCountFile);
	if($header) {
		chomp(my $line = <$reader>);
		print "$line\n";
	}
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		next if($line =~ /^__/);
		my ($gene, @countList) = split(/\t/, $line, -1);
		push(@geneCountListList, [$gene, @countList]);
		$totalCountList[$_] += $countList[$_] foreach(0 .. $#countList);
	}
	close($reader);
}
foreach(@geneCountListList) {
	my ($gene, @countList) = @$_;
	my @cpmList = map {$countList[$_] / ($totalCountList[$_] / 1000000)} 0 .. $#countList;
	print join("\t", $gene, @cpmList), "\n";
}
