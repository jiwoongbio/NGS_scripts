#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	's' => \(my $sorted = ''),
	'H=i' => \(my $headAddLength = 0),
	'T=i' => \(my $tailAddLength = 0),
	'v=s' => \(my $variantFile = ''),
);
my ($regionFile, $fastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($fastaFile);
my @nameList = ();
my %nameRegionListHash = ();
open(my $reader, $sorted ? $regionFile : "sort -t '\t' -k5 -k1,1 -k2,2n -k3,3n -k4,4 $regionFile |");
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome, $start, $end, $strand, $name) = split(/\t/, $line, 5);
	push(@nameList, $name) unless(defined($nameRegionListHash{$name}));
	push(@{$nameRegionListHash{$name}}, [$chromosome, $start, $end, $strand]);
}
close($reader);
foreach my $name (@nameList) {
	my @regionList = @{$nameRegionListHash{$name}};
	if($headAddLength > 0) {
		$regionList[0]->[1] -= $headAddLength if($regionList[0]->[3] eq '+');
		$regionList[-1]->[2] += $headAddLength if($regionList[-1]->[3] eq '-');
	}
	if($tailAddLength > 0) {
		$regionList[0]->[1] -= $tailAddLength if($regionList[0]->[3] eq '-');
		$regionList[-1]->[2] += $tailAddLength if($regionList[-1]->[3] eq '+');
	}
	my $sequence = '';
	foreach(@regionList) {
		my ($chromosome, $start, $end, $strand) = @$_;
		if($strand eq '+') {
			$sequence = $sequence . getSequence($chromosome, $start, $end);
		}
		if($strand eq '-') {
			$sequence = getReverseComplementarySequence(getSequence($chromosome, $start, $end)) . $sequence;
		}
	}
	print ">$name\n";
	print "$sequence\n";
}

sub getSequence {
	my ($chromosome, $start, $end) = @_;
	my $sequence = uc($db->seq($chromosome, $start, $end));
	if($variantFile ne '') {
		my @variantList = ();
		open(my $reader, "tabix $variantFile $chromosome:$start-$end |");
		while(my $line = <$reader>) {
			chomp($line);
			my ($chromosome, $position, $refBase, $altBase) = split(/\t/, $line, -1);
			($chromosome, $position, $refBase, $altBase) = extendIndel($chromosome, $position, $refBase, $altBase);
			next if($start > $position);
			next if($start > $position + length($refBase) - 1);
			next if($end < $position);
			next if($end < $position + length($refBase) - 1);
			my ($startIndex, $endIndex) = ($position - $start, $position + length($refBase) - 1 - $start);
			push(@variantList, [$startIndex, $endIndex, $chromosome, $position, $refBase, $altBase]);
		}
		close($reader);
		@variantList = sort {$b->[0] <=> $a->[0] || $b->[1] <=> $a->[1]} @variantList;
		my @baseList = split(//, $sequence);
		foreach my $index (0 .. $#variantList) {
			my ($startIndex, $endIndex, $chromosome, $position, $refBase, $altBase) = @{$variantList[$index]};
			if(grep {($startIndex == $variantList[$_]->[0] && $endIndex == $variantList[$_]->[1]) || ($startIndex <= $variantList[$_]->[1] && $variantList[$_]->[0] <= $endIndex)} grep {$_ != $index} 0 .. $#variantList) {
				print STDERR join("\t", 'overlap variant', $chromosome, $position, $refBase, $altBase), "\n";
			} elsif(join('', map {$baseList[$_]} $startIndex .. $endIndex) ne $refBase) {
				print STDERR join("\t", 'error variant', $chromosome, $position, $refBase, $altBase), "\n";
			} else {
				$baseList[$_] = '' foreach($startIndex .. $endIndex);
				$baseList[$startIndex] = $altBase . $baseList[$startIndex];
				print STDERR join("\t", 'apply variant', $chromosome, $position, $refBase, $altBase), "\n";
			}
		}
		$sequence = join('', @baseList);
	}
	return $sequence;
}

sub extendIndel {
	my ($chromosome, $position, $refBase, $altBase) = @_;
	if($refBase ne $altBase) {
		while($refBase =~ /^$altBase/ || $altBase =~ /^$refBase/) {
			my $extBase = uc($db->seq($chromosome, ($_ = $position + length($refBase)), $_));
			$refBase = "$refBase$extBase";
			$altBase = "$altBase$extBase";
		}
		while(substr($refBase, -1, 1) eq substr($altBase, -1, 1) && !($refBase =~ /^$altBase/ || $altBase =~ /^$refBase/)) {
			substr($refBase, -1, 1, '');
			substr($altBase, -1, 1, '');
		}
		while($refBase =~ /$altBase$/ || $altBase =~ /$refBase$/) {
			my $extBase = uc($db->seq($chromosome, ($position = $position - 1), $position));
			$refBase = "$extBase$refBase";
			$altBase = "$extBase$altBase";
		}
		while(substr($refBase, 0, 1) eq substr($altBase, 0, 1) && !($refBase =~ /$altBase$/ || $altBase =~ /$refBase$/)) {
			$position = $position + 1;
			substr($refBase, 0, 1, '');
			substr($altBase, 0, 1, '');
		}
	}
	return ($chromosome, $position, $refBase, $altBase);
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
