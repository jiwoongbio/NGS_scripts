#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'l=i' => \(my $lineLength = 80),
	'b' => \(my $binary = ''),
	'r' => \(my $includeReference = ''),
	'a=s' => \(my $ancestorSample = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl snv.table.pl [options] cohort.vcf

Options: -h       display this help message
         -l INT   line length [$lineLength]
         -b       binary
         -r       include reference
         -a       ancestor sample

EOF
}
my ($vcfFile) = @ARGV;
my @sampleList = ();
my @sequenceList = ();
my $ancestorSampleIndex = '';
{
	open(my $reader, ($vcfFile =~ /\.gz$/ ? "gzip -dc $vcfFile |" : $vcfFile));
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO//) {
			@sampleList = split(/\t/, $line, -1) if($line =~ s/^\tFORMAT\t//);
			if($ancestorSample ne '') {
				my %sampleIndexHash = map {$sampleList[$_] => $_} 0 .. $#sampleList;
				$ancestorSampleIndex = $sampleIndexHash{$ancestorSample};
			}
			next;
		}
		next if($line =~ /^#/);
		my %tokenHash = ();
		(@tokenHash{'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}, my @genotypeList) = map {$_ eq '.' ? '' : $_} split(/\t/, $line, -1);
		next if($tokenHash{'FILTER'} ne 'PASS');
		if(defined($tokenHash{'FORMAT'})) {
			my @alleleHashList = ();
			my @formatList = split(/:/, $tokenHash{'FORMAT'}, -1);
			foreach my $index (0 .. $#genotypeList) {
				my %tokenHash = ();
				@tokenHash{@formatList} = split(/:/, $genotypeList[$index], -1);
				$alleleHashList[$index]->{$_} = 1 foreach(grep {$_ ne '.'} split(/[\/|]/, $tokenHash{'GT'}));
			}
			if($binary eq '') {
				next if(grep {scalar(keys %$_) != 1} @alleleHashList);
				my @alleleList = map {keys %$_} @alleleHashList;
				@alleleList = (@alleleList, 0) if($includeReference);
				next if(scalar(unique(@alleleList)) == 1);
				my @alleleBaseList = ($tokenHash{'REF'}, $tokenHash{'ALT'} eq '' ? ('') : split(/,/, $tokenHash{'ALT'}, -1));
				@alleleBaseList = @alleleBaseList[@alleleList];
				next if(grep {$alleleBaseList[$_] !~ /^[ACGT]$/} 0 .. $#alleleList);
				$sequenceList[$_] .= $alleleBaseList[$_] foreach(0 .. $#alleleList);
			} else {
				next if(grep {scalar(keys %$_) == 0} @alleleHashList);
				foreach my $allele (unique(grep {$_ != 0} map {keys %$_} @alleleHashList)) {
					my @binaryList = map {$alleleHashList[$_]->{$allele} ? 1 : 0} 0 .. $#genotypeList;
					next if($ancestorSampleIndex ne '' && $binaryList[$ancestorSampleIndex] != 0);
					@binaryList = (@binaryList, 0) if($includeReference);
					next if(scalar(unique(@binaryList)) == 1);
					$sequenceList[$_] .= $binaryList[$_] foreach(0 .. $#binaryList);
				}
			}
		}
	}
	close($reader);
}
@sampleList = (@sampleList, 'reference') if($includeReference);
if(@sequenceList) {
	foreach my $index (0 .. $#sampleList) {
		print ">$sampleList[$index]\n";
		printSequence($sequenceList[$index]);
	}
}

sub printSequence {
	my ($sequence) = @_;
	for(my $index = 0; $index < length($sequence); $index += $lineLength) {
		print substr($sequence, $index, $lineLength), "\n";
	}
}

sub unique {
	my @sortedTokenList = sort @_;
	@sortedTokenList = @sortedTokenList[0, grep {$sortedTokenList[$_ - 1] ne $sortedTokenList[$_]} 1 .. $#sortedTokenList] if(@sortedTokenList);
	return @sortedTokenList;
}
