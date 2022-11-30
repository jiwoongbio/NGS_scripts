#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($k, @fastaFileList) = @ARGV;
foreach my $fastaFile (@fastaFileList) {
	open(my $reader, ($fastaFile =~ /\.gz$/ ? "gzip -dc $fastaFile |" : $fastaFile));
	my ($name, $sequence) = ('', '');
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			foreach my $index (0 .. length($sequence) - $k) {
				my $kmer = substr($sequence, $index, $k);
				print join("\t", $name, $index, $kmer), "\n";
			}
			$name = $1;
			$sequence = '';
		} else {
			$sequence .= $line;
		}
	}
	{
		foreach my $index (0 .. length($sequence) - $k) {
			my $kmer = substr($sequence, $index, $k);
			print join("\t", $name, $index, $kmer), "\n";
		}
		$name = $1;
		$sequence = '';
	}
	close($reader);
}
