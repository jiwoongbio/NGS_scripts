#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'r=s' => \(my $referenceFastaFile = ''),
	'a=s' => \(my $abundanceFile = ''),
);
my ($k, @fastaFileList) = @ARGV;
my %kmerCountHash = ();
if($referenceFastaFile ne '') {
	open(my $reader, ($referenceFastaFile =~ /\.gz$/ ? "gzip -dc $referenceFastaFile |" : $referenceFastaFile));
	my ($name, $sequence) = ('', '');
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			foreach my $index (0 .. length($sequence) - $k) {
				my $kmer = substr($sequence, $index, $k);
				$kmerCountHash{$kmer} = 0;
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
			$kmerCountHash{$kmer} = 0;
		}
		$name = '';
		$sequence = '';
	}
	close($reader);
}
my %abundanceHash = ();
if($abundanceFile ne '') {
	open(my $reader, ($abundanceFile =~ /\.gz$/ ? "gzip -dc $abundanceFile |" : $abundanceFile));
	while(my $line = <$reader>) {
		chomp($line);
		my ($name, $abundance) = split(/\t/, $line, -1);
		$abundanceHash{$name} = $abundance;
	}
	close($reader);
}
foreach my $fastaFile (@fastaFileList) {
	open(my $reader, ($fastaFile =~ /\.gz$/ ? "gzip -dc $fastaFile |" : $fastaFile));
	my ($name, $sequence) = ('', '');
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			my $abundance = 1;
			$abundance = $abundanceHash{$name} if($abundanceFile ne '');
			foreach my $index (0 .. length($sequence) - $k) {
				my $kmer = substr($sequence, $index, $k);
				if($referenceFastaFile ne '') {
					$kmerCountHash{$kmer} += $abundance if(defined($kmerCountHash{$kmer}));
				} else {
					$kmerCountHash{$kmer} += $abundance;
				}
			}
			$name = $1;
			$sequence = '';
		} else {
			$sequence .= $line;
		}
	}
	{
		my $abundance = 1;
		$abundance = $abundanceHash{$name} if($abundanceFile ne '');
		foreach my $index (0 .. length($sequence) - $k) {
			my $kmer = substr($sequence, $index, $k);
			if($referenceFastaFile ne '') {
				$kmerCountHash{$kmer} += $abundance if(defined($kmerCountHash{$kmer}));
			} else {
				$kmerCountHash{$kmer} += $abundance;
			}
		}
		$name = $1;
		$sequence = '';
	}
	close($reader);
}
foreach my $kmer (sort keys %kmerCountHash) {
	print join("\t", $kmer, $kmerCountHash{$kmer}), "\n" if($kmerCountHash{$kmer});
}
