#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

my @sampleList = ();
GetOptions(
	's=s' => \@sampleList,
	'S=s' => \(my $sampleFile = 'sample.txt'),
	'h' => \(my $header = ''),
);
@sampleList = map {split(/,/, $_)} @sampleList;

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
	['table.variant_count.not_dbSNP_not_assembly_variant.txt', 'Not in dbSNP and not assembly variant'],
	['mutect.table.variant_count.txt', 'Somatic mutation calling'],
	['mutect.table.variant_count.not_Common.txt', 'Not in dbSNP Common subset'],
	['mutect.table.variant_count.not_dbSNP.txt', 'Not in dbSNP'],
	['mutect.table.variant_count.not_dbSNP_not_assembly_variant.txt', 'Not in dbSNP and not assembly variant'],
);
@suffixNameList = ((map {["remocon.$_->[0]", $_->[1]]} @suffixNameList), @suffixNameList);
my @tandemRepeatSuffixList = ('remocon.mutect.tandemRepeat.count.txt', 'remocon.tandemRepeat.count.txt', 'mutect.tandemRepeat.count.txt', 'tandemRepeat.count.txt');

if(my @variantFileList = @ARGV) {
	my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2 -k3,6 | uniq | cut -f1,2 | uniq -c");
	foreach my $variantFile (@variantFileList) {
		my %columnIndexHash = ();
		open(my $reader, $variantFile);
		if($header ne '') {
			chomp(my $line = <$reader>);
			$line =~ s/^#//;
			my @columnList = split(/\t/, $line, -1);
			%columnIndexHash = getColumnIndexHash(@columnList);
			$columnIndexHash{'chromosome'} = $columnIndexHash{'Chromosome'} unless(defined($columnIndexHash{'chromosome'}));
			$columnIndexHash{'position'} = $columnIndexHash{'Position'} unless(defined($columnIndexHash{'position'}));
			$columnIndexHash{'refBase'} = $columnIndexHash{'Reference base'} unless(defined($columnIndexHash{'refBase'}));
			$columnIndexHash{'altBase'} = $columnIndexHash{'Alternate base'} unless(defined($columnIndexHash{'altBase'}));
			$columnIndexHash{'sample'} = $columnIndexHash{'Sample'} unless(defined($columnIndexHash{'sample'}));
			$columnIndexHash{'genotype'} = $columnIndexHash{'Genotype'} unless(defined($columnIndexHash{'genotype'}));
			$columnIndexHash{'genotype'} = $columnIndexHash{'Tumor genotype'} unless(defined($columnIndexHash{'genotype'}));
			$columnIndexHash{'gene'} = $columnIndexHash{'Gene name'} unless(defined($columnIndexHash{'gene'}));
			$columnIndexHash{'region'} = $columnIndexHash{'Region type'} unless(defined($columnIndexHash{'region'}));
			$columnIndexHash{'mutation'} = $columnIndexHash{'Mutation class'} unless(defined($columnIndexHash{'mutation'}));
		} else {
			my @columnList = ('chromosome', 'position', 'refBase', 'altBase', 'sample', 'genotype', 'gene', 'region', 'mutation');
			%columnIndexHash = getColumnIndexHash(@columnList);
		}
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^#/);
			my @tokenList = split(/\t/, $line, -1);
			my ($chromosome, $position, $refBase, $altBase, $sample, $genotype, $gene, $region, $mutation) = map {defined($_) ? $tokenList[$_] : $_} @columnIndexHash{'chromosome', 'position', 'refBase', 'altBase', 'sample', 'genotype', 'gene', 'region', 'mutation'};
			$sample = $sampleList[0] if(!defined($sample) && scalar(@sampleList) == 1);
			foreach my $column (@columnList) {
				my $expression = $columnExpressionHash{$column};
				print $writer join("\t", $sample, $column, $chromosome, $position, $refBase, $altBase), "\n" if(eval($expression));
			}
		}
		close($reader);
	}
	close($writer);
	my %sampleColumnCountHash = ();
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^ *([0-9]+) (.*)$/) {
			my $count = $1;
			my ($sample, $column) = split(/\t/, $2, -1);
			$sampleColumnCountHash{$sample}->{$column} = $count;
		}
	}
	close($reader);
	waitpid($pid, 0);
	@sampleList = sort keys %sampleColumnCountHash unless(@sampleList);
	print join("\t", 'sample', @sampleList), "\n";
	foreach my $column (@columnList) {
		print join("\t", $column, map {defined($_) ? $_ : 0} map {$_->{$column}} @sampleColumnCountHash{@sampleList}), "\n";
	}
} elsif(-s $sampleFile) {
	chomp(@sampleList = `cat $sampleFile`) unless(@sampleList);
	my @sampleDirectoryPrefixList = ();
	foreach my $sample (@sampleList) {
		chomp(my $directory = `readlink -f $sample`);
		(my $prefix = $directory) =~ s/^.*\///;
		push(@sampleDirectoryPrefixList, [$sample, "$directory/$prefix"]);
	}
	my @availableSuffixNameList = ();
	foreach(@suffixNameList) {
		my ($suffix, $name) = @$_;
		push(@availableSuffixNameList, [$suffix, $name]) if(grep {-s "$_->[1].$suffix"} @sampleDirectoryPrefixList);
	}
	my @availableTandemRepeatSuffixList = ();
	foreach my $suffix (@tandemRepeatSuffixList) {
		push(@availableTandemRepeatSuffixList, $suffix) if(grep {-s "$_->[1].$suffix"} @sampleDirectoryPrefixList);
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
	foreach(@sampleDirectoryPrefixList) {
		my ($sample, $directoryPrefix) = @$_;
		my @tokenList = ($sample);
		foreach(@availableSuffixNameList) {
			my ($suffix, $name) = @$_;
			my %tokenHash = ();
			if(-s "$directoryPrefix.$suffix") {
				open(my $reader, "$directoryPrefix.$suffix");
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
			if(-s (my $file = "$directoryPrefix.$availableTandemRepeatSuffixList[0]")) {
				chomp(my $count = `cat $file`);
				push(@tokenList, $count);
			} else {
				push(@tokenList, '');
			}
		}
		print join("\t", @tokenList), "\n";
	}
}

sub getColumnIndexHash {
	my (@columnList) = @_;
	my %columnIndexHash = ();
	@columnIndexHash{@columnList} = 0 .. $#columnList;
	foreach my $column ('chromosome', 'Chromosome', 'CHROM', 'chr') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'chromosome'} = $index;
			last;
		}
	}
	foreach my $column ('position', 'Position', 'POS') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'position'} = $index;
			last;
		}
	}
	foreach my $column ('refBase', 'haplotypeReference', 'Reference base', 'REF') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'refBase'} = $index;
			last;
		}
	}
	foreach my $column ('altBase', 'haplotypeAlternate', 'Alternate base', 'ALT') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'altBase'} = $index;
			last;
		}
	}
	foreach my $column ('sample', 'Sample') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'sample'} = $index;
			last;
		}
	}
	foreach my $column ('genotype', 'Genotype') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'genotype'} = $index;
			last;
		}
	}
	foreach my $column ('gene', 'Gene name') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'gene'} = $index;
			last;
		}
	}
	foreach my $column ('region', 'Region type') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'region'} = $index;
			last;
		}
	}
	foreach my $column ('mutation', 'Mutation class') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'mutation'} = $index;
			last;
		}
	}
	return %columnIndexHash;
}

sub getHaplotypeList {
	my ($genotype) = @_;
	my %haplotypeHash = ();
	$haplotypeHash{$_} = 1 foreach(split(/[\/|]/, $genotype));
	my @haplotypeList = sort {$a <=> $b} grep {$_ ne '.'} keys %haplotypeHash;
	return @haplotypeList;
}
