#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	's=s' => \(my $samples = ''),
);

my %columnExpressionHash = ();
$columnExpressionHash{'total'} = '1';

$columnExpressionHash{'homozygous'}   = 'scalar(getHaplotypeList($genotype)) == 1';
$columnExpressionHash{'heterozygous'} = 'scalar(getHaplotypeList($genotype)) > 1';

$columnExpressionHash{'transition'}   = 'grep {$_->[0] eq $refBase && $_->[1] eq $altBase} ["A", "G"], ["C", "T"], ["G", "A"], ["T", "C"]';
$columnExpressionHash{'transversion'} = 'grep {$_->[0] eq $refBase && $_->[1] eq $altBase} ["A", "C"], ["A", "T"], ["C", "A"], ["C", "G"], ["G", "C"], ["G", "T"], ["T", "A"], ["T", "G"]';

my @regionList = ('CDS', 'utr5', 'utr3', 'intron', 'non_coding_exon', 'non_coding_intron');
$columnExpressionHash{$_} = "\$region =~ /$_/" foreach(@regionList);
$columnExpressionHash{'exon'} = join(' || ', map {$columnExpressionHash{$_}} ('CDS', 'utr5', 'utr3', 'non_coding_exon'));

my @mutationList = ('silent', 'missense', 'inframe', 'frameshift', 'nonsense', 'readthrough', 'startcodon', 'splicing', 'junction');
$columnExpressionHash{$_} = "\$mutation =~ /$_/" foreach(@mutationList);
$columnExpressionHash{'NS/SS/I'} = join(' || ', map {$columnExpressionHash{$_}} grep {$_ ne 'silent'} @mutationList);

my @columnList = ('total', 'homozygous', 'heterozygous', 'transition', 'transversion', 'exon', @regionList, @mutationList, 'NS/SS/I');

my @suffixNameList = (
	['table.variant_count.txt', 'Variant calling'],
	['table.variant_count.not_Common.txt', 'Not in dbSNP Common subset'],
	['table.variant_count.not_dbSNP.txt', 'Not in dbSNP'],
	['mutect.table.variant_count.txt', 'Somatic mutation calling'],
	['mutect.table.variant_count.not_Common.txt', 'Not in dbSNP Common subset'],
	['mutect.table.variant_count.not_dbSNP.txt', 'Not in dbSNP'],
);
@suffixNameList = ((map {["remocon.$_->[0]", $_->[1]]} @suffixNameList), @suffixNameList);
my @tandemRepeatSuffixList = ('remocon.mutect.tandemRepeat.count.txt', 'remocon.tandemRepeat.count.txt', 'mutect.tandemRepeat.count.txt', 'tandemRepeat.count.txt');

if(my @variantFileList = @ARGV) {
my @sampleList = ();
	if($samples eq '') {
		my %sampleHash = ();
		foreach my $variantFile (@variantFileList) {
			open(my $reader, $variantFile);
			while(my $line = <$reader>) {
				chomp($line);
				next if($line =~ /^#/);
				my ($chromosome, $position, $refBase, $altBase, $sample, $genotype, $gene, $region, $mutation) = split(/\t/, $line, -1);
				$sampleHash{$sample} = 1;
			}
			close($reader);
		}
		@sampleList = sort keys %sampleHash;
	} else {
		@sampleList = split(/,/, $samples);
	}
	print join("\t", 'sample', @sampleList), "\n";
	foreach my $column (@columnList) {
		my $expression = $columnExpressionHash{$column};
		my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,5 | uniq | cut -f1 | uniq -c");
		foreach my $variantFile (@variantFileList) {
			open(my $reader, $variantFile);
			while(my $line = <$reader>) {
				chomp($line);
				next if($line =~ /^#/);
				my ($chromosome, $position, $refBase, $altBase, $sample, $genotype, $gene, $region, $mutation) = split(/\t/, $line, -1);
				print $writer join("\t", $sample, $chromosome, $position, $refBase, $altBase), "\n" if(eval($expression));
			}
			close($reader);
		}
		close($writer);
		my %sampleCountHash = ();
		while(my $line = <$reader>) {
			chomp($line);
			$sampleCountHash{$2} = $1 if($line =~ /^ *([0-9]+) (.*)$/);
		}
		close($reader);
		waitpid($pid, 0);
		print join("\t", $column, map {defined($_) ? $_ : 0} @sampleCountHash{@sampleList}), "\n";
	}
} elsif(-s 'sample.txt') {
	chomp(my @sampleList = `cat sample.txt`);
	my @availableSuffixNameList = ();
	foreach(@suffixNameList) {
		my ($suffix, $name) = @$_;
		push(@availableSuffixNameList, [$suffix, $name]) if(grep {-s "$_/$_.$suffix"} @sampleList);
	}
	my @availableTandemRepeatSuffixList = ();
	foreach my $suffix (@tandemRepeatSuffixList) {
		push(@availableTandemRepeatSuffixList, $suffix) if(grep {-s "$_/$_.$suffix"} @sampleList);
	}
	if(scalar(@availableSuffixNameList) > 1) {
		my @tokenList = ('');
		foreach(@availableSuffixNameList) {
			my ($suffix, $name) = @$_;
			push(@tokenList, $name, ('') x (scalar(@columnList) - 1));
		}
		push(@tokenList, '') if(@availableTandemRepeatSuffixList);
		print join("\t", @tokenList), "\n";
	}
	{
		my @tokenList = ('sample');
		foreach(@availableSuffixNameList) {
			my ($suffix, $name) = @$_;
			push(@tokenList, @columnList);
		}
		push(@tokenList, 'tandem repeat') if(@availableTandemRepeatSuffixList);
		print join("\t", @tokenList), "\n";
	}
	foreach my $sample (@sampleList) {
		my @tokenList = ($sample);
		foreach(@availableSuffixNameList) {
			my ($suffix, $name) = @$_;
			my %tokenHash = ();
			if(-s "$sample/$sample.$suffix") {
				open(my $reader, "$sample/$sample.$suffix");
				while(my $line = <$reader>) {
					chomp($line);
					my @tokenList = split(/\t/, $line, -1);
					$tokenHash{$tokenList[0]} = $tokenList[1];
				}
				close($reader);
			}
			push(@tokenList, map {defined($_) ? $_ : ''} @tokenHash{@columnList});
		}
		if(@availableTandemRepeatSuffixList) {
			if(-s (my $file = "$sample/$sample.$availableTandemRepeatSuffixList[0]")) {
				chomp(my $count = `cat $file`);
				push(@tokenList, $count);
			} else {
				push(@tokenList, '');
			}
		}
		print join("\t", @tokenList), "\n";
	}
}

sub getHaplotypeList {
	my ($genotype) = @_;
	my %haplotypeHash = ();
	$haplotypeHash{$_} = 1 foreach(split(/[\/|]/, $genotype));
	my @haplotypeList = sort {$a <=> $b} grep {$_ ne '.'} keys %haplotypeHash;
	return @haplotypeList;
}
