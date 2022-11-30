#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($pileupFile, $referenceFastaFile, $variantFile) = @ARGV;
my %variantBaseHash = ();
if(defined($variantFile)) {
	my @chromosomeList = getChromosomeList();
	my %chromosomeIndexHash = ();
	@chromosomeIndexHash{@chromosomeList} = 0 .. scalar(@chromosomeList) - 1;
	open(my $reader, ($variantFile =~ /\.gz$/ ? "gzip -dc $variantFile |" : $variantFile));
	sub getSingleNucleotideVariant {
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^#/);
			my @tokenList = split(/\t/, $line, -1);
			return \@tokenList if($tokenList[2] =~ /^[ACGT]$/ && $tokenList[3] =~ /^[ACGT]$/);
		}
		return '';
	}
	my @variantList = ();
	my $variant = getSingleNucleotideVariant();
	sub setVariantChromosomePosition {
		my ($chromosome, $position) = @_;
		return if(@variantList && $variantList[0]->[0] eq $chromosome && $variantList[0]->[1] == $position);
		@variantList = ();
		while($variant && $chromosomeIndexHash{$variant->[0]} < $chromosomeIndexHash{$chromosome}) {
			$variant = getSingleNucleotideVariant();
		}
		while($variant && $variant->[0] eq $chromosome && $variant->[1] < $position) {
			$variant = getSingleNucleotideVariant();
		}
		while($variant && $variant->[0] eq $chromosome && $variant->[1] == $position) {
			push(@variantList, [$variant]);
			$variant = getSingleNucleotideVariant();
		}
		$variantBaseHash{$_->[3]} = 1 foreach(@variantList);
	}
	sub closeVariantFileReader {
		close($reader) if($reader);
	}
}
my %baseBaseCountHash = ();
open(my $reader, ($pileupFile =~ /\.gz$/ ? "gzip -dc $pileupFile |" : $pileupFile));
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split("\t", $line, -1);
	my ($chromosome, $position, $refBase) = @tokenList;
	setVariantChromosomePosition($chromosome, $position);
	$refBase =~ tr/a-z/A-Z/;
	my %readBaseCountHash = ();
	for(my $index = 3; $index < scalar(@tokenList); $index += 3) {
		my ($readDepth, $readBases, $baseQualities) = @tokenList[$index .. $index + 2];
		foreach my $readBase (getReadBaseList($readBases)) {
			$readBase =~ s/^\^.//;
			$readBase =~ s/\$$//;
			next if($readBase eq '>');
			next if($readBase eq '<');
			next if($readBase =~ /^\*/);
			$readBase =~ s/^[.,]/$refBase/;
			$readBase =~ tr/a-z/A-Z/;
			$readBaseCountHash{$readBase} += 1 if($readBase =~ /^[ACGT]$/);
		}
	}
	$baseBaseCountHash{$refBase}->{$_} += $readBaseCountHash{$_} foreach(grep {!$variantBaseHash{$_}} keys %readBaseCountHash);
}
close($reader);
foreach my $refBase ('A', 'C', 'G', 'T') {
	foreach my $readBase ('A', 'C', 'G', 'T') {
		if(my $count = $baseBaseCountHash{$refBase}->{$readBase}) {
			print join("\t", $refBase, $readBase, $count), "\n";
		} else {
			print join("\t", $refBase, $readBase, 0), "\n";
		}
	}
}

sub getReadBaseList {
	my ($readBases) = @_;
	my @readBaseList = ();
	while($readBases ne '') {
		my $readBase = '';
		$readBase .= $1 if($readBases =~ s/^(\^.)//);
		$readBases =~ s/^([.A-Z>,a-z<*])//;
		$readBase .= $1;
		while($readBases =~ s/^([+-])([0-9]+)//) {
			$readBase .= "$1$2";
			$readBase .= substr($readBases, 0, $2, '');
		}
		$readBase .= $1 if($readBases =~ s/^(\$)//);
		push(@readBaseList, $readBase);
	}
	return @readBaseList;
}

sub getChromosomeList {
	my @chromosomeList = ();
	if(my $faiFile = `find $referenceFastaFile.fai -newer $referenceFastaFile 2> /dev/null`) {
		chomp($faiFile);
		open(my $reader, $faiFile);
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line);
			push(@chromosomeList, $tokenList[0]);
		}
		close($reader);
	} else {
		open(my $reader, $referenceFastaFile);
		while(my $line = <$reader>) {
			chomp($line);
			push(@chromosomeList, $1) if($line =~ /^>(\S*)/);
		}
		close($reader);
	}
	return @chromosomeList;
}
