#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;
use List::Util qw(sum);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'd=i' => \(my $minimumReadBaseDepth = 0),
);
my ($pileupFile, $referenceFastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
open(my $reader, ($pileupFile =~ /\.gz$/ ? "gzip -dc $pileupFile |" : $pileupFile));
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split("\t", $line, -1);
	my ($chromosome, $position, $refBase) = @tokenList;
	$refBase =~ tr/a-z/A-Z/;
	my @readDepthList = ();
	my %readBaseDepthListHash = ();
	for(my $index = 3; $index < scalar(@tokenList); $index += 3) {
		my ($readDepth, $readBases, $baseQualities) = @tokenList[$index .. $index + 2];
		push(@readDepthList, $readDepth);
		foreach my $readBase (getReadBaseList($readBases)) {
			$readBase =~ s/^\^.//;
			$readBase =~ s/\$$//;
			$readBase =~ s/^[.,]/$refBase/;
			$readBase =~ tr/a-z/A-Z/;
			$readBaseDepthListHash{$readBase}->[($index / 3) - 1] += 1;
		}
	}
	foreach my $readBase (grep {$_ ne $refBase && $_ ne '>' && $_ ne '<' && $_ ne '*'} sort keys %readBaseDepthListHash) {
		if(sum(0, grep {defined} @{$readBaseDepthListHash{$readBase}}) >= $minimumReadBaseDepth) {
			printVariant($chromosome, $position, $refBase, $readBase, map {defined($_) ? $_ : 0} map {($readDepthList[$_], $readBaseDepthListHash{$readBase}->[$_])} 0 .. $#readDepthList);
		}
	}
}
close($reader);

sub printVariant {
	my ($chromosome, $position, $refBase, $altBase, @depthList) = @_;
	if($altBase =~ s/([+-])([0-9]+)([ACGTNacgtn]+)$//) {
		$refBase .= $3 if($1 eq '-');
		$altBase .= $3 if($1 eq '+');
	}
	print join("\t", extendIndel($chromosome, $position, $refBase, $altBase), map {($depthList[$_ * 2], $depthList[$_ * 2 + 1], $depthList[$_ * 2] == 0 ? 0 : $depthList[$_ * 2 + 1] / $depthList[$_ * 2])} 0 .. scalar(@depthList) / 2 - 1), "\n";
}

sub getReadBaseList {
	my ($readBases) = @_;
	my @readBaseList = ();
	while($readBases ne '') {
		my $readBase = '';
		$readBase .= $1 if($readBases =~ s/^(\^.)//);
		$readBases =~ s/^([.ACGTN>,acgtn<*])//;
		$readBase .= $1;
		if($readBases =~ s/^([+-])([0-9]+)//) {
			$readBase .= "$1$2";
			$readBases =~ s/^([ACGTNacgtn]{$2})//;
			$readBase .= $1;
		}
		$readBase .= $1 if($readBases =~ s/^(\$)//);
		push(@readBaseList, $readBase);
	}
	return @readBaseList;
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
