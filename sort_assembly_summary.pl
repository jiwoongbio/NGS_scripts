#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;

my @refseqCategoryList = ('reference genome', 'representative genome', 'na');
my %refseqCategoryIndexHash = map {$refseqCategoryList[$_] => $_} 0 .. $#refseqCategoryList;

my @assemblyLevelList = ('Complete Genome', 'Chromosome', 'Scaffold', 'Contig');
my %assemblyLevelIndexHash = map {$assemblyLevelList[$_] => $_} 0 .. $#assemblyLevelList;

my (@assemblySummaryFileList) = @ARGV;
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1n -k2,2n -k3,3n -k4,4n | cut -f5-");
my @outputColumnList = ();
{
	foreach my $assemblySummaryFileIndex (0 .. $#assemblySummaryFileList) {
		my $assemblySummaryFile = $assemblySummaryFileList[$assemblySummaryFileIndex];
		open(my $reader, ($assemblySummaryFile =~ /\.gz$/ ? "gzip -dc $assemblySummaryFile |" : $assemblySummaryFile)) or die "Can't open '$assemblySummaryFile': $!";
		my $line;
		chomp($line = <$reader>);
		print $line, "\n" if($assemblySummaryFileIndex == 0);
		chomp($line = <$reader>);
		print $line, "\n" if($assemblySummaryFileIndex == 0);
		$line =~ s/^# //;
		my @columnList = split(/\t/, $line, -1);
		@outputColumnList = @columnList if($assemblySummaryFileIndex == 0);
		while($line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			my $refseqCategoryIndex = $refseqCategoryIndexHash{$tokenHash{'refseq_category'}};
			my $assemblyLevelIndex = $assemblyLevelIndexHash{$tokenHash{'assembly_level'}};
			(my $assemblyAccessionNumber = $tokenHash{'assembly_accession'}) =~ s/^[^0-9]*//;
			print $writer join("\t", $refseqCategoryIndex, $assemblyLevelIndex, $assemblySummaryFileIndex, $assemblyAccessionNumber, @tokenHash{@outputColumnList}), "\n";
		}
		close($reader);
	}
}
close($writer);
{
	my %assemblyAccessionHash = ();
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@outputColumnList} = split(/\t/, $line, -1);
		print $line, "\n" unless($assemblyAccessionHash{$tokenHash{'assembly_accession'}});
		$assemblyAccessionHash{$tokenHash{'assembly_accession'}} = 1;
		$assemblyAccessionHash{$tokenHash{'gbrs_paired_asm'}} = 1 if($tokenHash{'paired_asm_comp'} eq 'identical');
	}
}
close($reader);
waitpid($pid, 0);
