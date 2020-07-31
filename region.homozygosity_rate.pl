#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;

my ($regionFile, $referenceFastaFile, $vcfFile, $sample) = @ARGV;
{
	my @chromosomeList = getChromosomeList();
	my %chromosomeIndexHash = ();
	@chromosomeIndexHash{@chromosomeList} = 0 .. $#chromosomeList;
	open(my $reader, ($vcfFile =~ /\.gz$/) ? "gzip -dc $vcfFile |" : $vcfFile);
	chomp(my $line = <$reader>);
	my @sampleList = ();
	my $tokenHash = getTokenHash();
	my @tokenHashList = ();
	sub getTokenHash {
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ s/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//) {
				@sampleList = split(/\t/, $line, -1);
				next;
			}
			next if($line =~ /^#/);
			my %tokenHash = ();
			(@tokenHash{'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}, my @genotypeList) = map {$_ eq '.' ? '' : $_} split(/\t/, $line, -1);
			next if($tokenHash{'FILTER'} ne 'PASS');
			if($tokenHash{'REF'} =~ /^[ACGT]$/ && $tokenHash{'ALT'} =~ /^[ACGT]$/) {
				my %sampleGenotypeHash = ();
				@sampleGenotypeHash{@sampleList} = @genotypeList;
				@tokenHash{split(/:/, $tokenHash{'FORMAT'}, -1)} = split(/:/, $sampleGenotypeHash{$sample}, -1);
				@tokenHash{'chromosome', 'start', 'end'} = @tokenHash{'CHROM', 'POS', 'POS'};
				my @alleleList = unique(grep {$_ ne '.'} split(/[\/|]/, $tokenHash{'GT'}));
				$tokenHash{'alleleCount'} = scalar(@alleleList);
				return \%tokenHash;
			}
		}
		return '';
	}
	sub getTokenHashList {
		my ($chromosome, $start, $end) = @_;
		@tokenHashList = grep {$_->{'chromosome'} eq $chromosome && $start <= $_->{'end'}} @tokenHashList;
		while($tokenHash && $chromosomeIndexHash{$tokenHash->{'chromosome'}} < $chromosomeIndexHash{$chromosome}) {
			$tokenHash = getTokenHash();
		}
		while($tokenHash && $tokenHash->{'chromosome'} eq $chromosome && $tokenHash->{'start'} <= $end) {
			push(@tokenHashList, $tokenHash) if($start <= $tokenHash->{'end'});
			$tokenHash = getTokenHash();
		}
		return grep {$_->{'start'} <= $end} @tokenHashList;
	}
	sub getHomozygosityRate {
		my ($chromosome, $start, $end) = @_;
		if(my @tokenHashList = getTokenHashList($chromosome, $start, $end)) {
			return scalar(grep {$_->{'alleleCount'} == 1} @tokenHashList) / scalar(@tokenHashList);
		} else {
			return '';
		}
	}
	sub closeVcfFileReader {
		close($reader);
	}
}
open(my $reader, ($regionFile =~ /\.gz$/) ? "gzip -dc $regionFile |" : $regionFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^#/) {
		print "$line\n";
	} else {
		my @tokenList = split(/\t/, $line, -1);
		my ($chromosome, $start, $end) = @tokenList;
		print join("\t", @tokenList, getHomozygosityRate($chromosome, $start, $end)), "\n";
	}
}
close($reader);
closeVcfFileReader();

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

sub unique {
	my @sortedTokenList = sort @_;
	return @sortedTokenList[0, grep {$sortedTokenList[$_ - 1] ne $sortedTokenList[$_]} 1 .. $#sortedTokenList];
}
