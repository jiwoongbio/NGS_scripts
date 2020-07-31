#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case);

my @codonList = ();
GetOptions(
	'h' => \(my $help = ''),
	'C=s' => \@codonList,
	'e' => \(my $extend = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   region.translate.pl [options] region.txt sequence.fasta > translation.fasta

Options: -h       display this help message
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 1 (standard)]
         -e       extend

EOF
}
{
	my %codonHash = (
		'TTT' => 'F', 'CTT' => 'L', 'ATT' => 'I', 'GTT' => 'V',
		'TTC' => 'F', 'CTC' => 'L', 'ATC' => 'I', 'GTC' => 'V',
		'TTA' => 'L', 'CTA' => 'L', 'ATA' => 'I', 'GTA' => 'V',
		'TTG' => 'L', 'CTG' => 'L', 'ATG' => 'M', 'GTG' => 'V',

		'TCT' => 'S', 'CCT' => 'P', 'ACT' => 'T', 'GCT' => 'A',
		'TCC' => 'S', 'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
		'TCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
		'TCG' => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',

		'TAT' => 'Y', 'CAT' => 'H', 'AAT' => 'N', 'GAT' => 'D',
		'TAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
		'TAA' => '*', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
		'TAG' => '*', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',

		'TGT' => 'C', 'CGT' => 'R', 'AGT' => 'S', 'GGT' => 'G',
		'TGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
		'TGA' => '*', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
		'TGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G',
	);
	$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} @codonList);

	sub translate {
		my ($sequence) = @_;
		return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. int(length($sequence) / 3) - 1);
	}
}
my ($regionFile, $fastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($fastaFile);
my %nameChromosomeStartEndStrandListHash = ();
open(my $reader, $regionFile);
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome, $start, $end, $strand, $name) = split(/\t/, $line, -1);
	push(@{$nameChromosomeStartEndStrandListHash{$name}}, [$chromosome, $start, $end, $strand]);
}
close($reader);
foreach my $name (sort keys %nameChromosomeStartEndStrandListHash) {
	my @chromosomeStartEndStrandList = @{$nameChromosomeStartEndStrandListHash{$name}};
	if($extend) {
		$chromosomeStartEndStrandList[0]->[1] = 1 if($chromosomeStartEndStrandList[0]->[3] eq '-');
		$chromosomeStartEndStrandList[-1]->[2] = $db->length($chromosomeStartEndStrandList[-1]->[0]) if($chromosomeStartEndStrandList[-1]->[3] eq '+');
	}
	next if(grep {$_->[1] eq '' || $_->[2] eq ''} @chromosomeStartEndStrandList);
	my $sequence = '';
	foreach(@chromosomeStartEndStrandList) {
		my ($chromosome, $start, $end, $strand) = @$_;
		if($strand eq '+') {
			$sequence = $sequence . uc($db->seq($chromosome, $start, $end));
		}
		if($strand eq '-') {
			$sequence = uc($db->seq($chromosome, $end, $start)) . $sequence;
		}
	}
	(my $translationSequence = translate($sequence)) =~ s/\*.*$//;
	if($translationSequence ne '') {
		print ">$name\n";
		print "$translationSequence\n";
	}
}
