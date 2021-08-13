#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'H' => \(my $header = ''),
	'g=s' => \(my $geneAttribute = 'gene_id'),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl RPKM.pl [options] gene.gtf gene.count.txt > RPKM.txt

Options: -h       display this help message
         -H       with header
         -g STR   gene attribute

EOF
}
my ($gtfFile, $geneCountFile) = @ARGV;
die "GTF file is not defined.\n" unless(defined($gtfFile));
die "$gtfFile is not readable.\n" unless(-r $gtfFile);
die "Gene count file is not defined.\n" unless(defined($geneCountFile));
die "$geneCountFile is not readable.\n" unless(-r $geneCountFile);
my %geneLengthHash = ();
{
	my @columnList = ($geneAttribute, 'chromosome', 'start', 'end');
	my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2 -k3,3n -k4,4n | uniq");
	{
		open(my $reader, $gtfFile);
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^#/);
			my %tokenHash = ();
			@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line, -1);
			$tokenHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g);
			print $writer join("\t", @tokenHash{@columnList}), "\n";
		}
		close($reader);
	}
	close($writer);
	my %geneChromosomeEndHash = ();
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		if(defined(my $end = $geneChromosomeEndHash{$tokenHash{$geneAttribute}}->{$tokenHash{'chromosome'}})) {
			if($end < $tokenHash{'end'}) {
				if($end < $tokenHash{'start'}) {
					$geneLengthHash{$tokenHash{$geneAttribute}} += $tokenHash{'end'} - $tokenHash{'start'} + 1;
				} else {
					$geneLengthHash{$tokenHash{$geneAttribute}} += $tokenHash{'end'} - $end;
				}
				$geneChromosomeEndHash{$tokenHash{$geneAttribute}}->{$tokenHash{'chromosome'}} = $tokenHash{'end'};
			}
		} else {
			$geneLengthHash{$tokenHash{$geneAttribute}} += $tokenHash{'end'} - $tokenHash{'start'} + 1;
			$geneChromosomeEndHash{$tokenHash{$geneAttribute}}->{$tokenHash{'chromosome'}} = $tokenHash{'end'};
		}
	}
	close($reader);
	waitpid($pid, 0);
}
my @geneCountListList = ();
my @totalCountList = ();
{
	open(my $reader, $geneCountFile);
	if($header) {
		chomp(my $line = <$reader>);
		print "$line\n";
	}
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		next if($line =~ /^__/);
		my ($gene, @countList) = split(/\t/, $line, -1);
		push(@geneCountListList, [$gene, @countList]);
		$totalCountList[$_] += $countList[$_] foreach(0 .. $#countList);
	}
	close($reader);
}
foreach(@geneCountListList) {
	my ($gene, @countList) = @$_;
	if($geneLengthHash{$gene}) {
		my @rpkmList = map {$countList[$_] / ($geneLengthHash{$gene} / 1000) / ($totalCountList[$_] / 1000000)} 0 .. $#countList;
		print join("\t", $gene, @rpkmList), "\n";
	} else {
		print STDERR "$gene not in GTF\n";
	}
}
