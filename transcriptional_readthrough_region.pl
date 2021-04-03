#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'g=s' => \(my $geneAttribute = 'gene_id'),
	't=s' => \(my $transcriptAttribute = 'transcript_id'),
);
my ($gtfFile, $referenceFastaFile, $maximumLength) = @ARGV;
my %chromosomeSequenceLengthHash = ();
{
	my $chromosome = '';
	open(my $reader, ($referenceFastaFile =~ /\.gz$/) ? "gzip -dc $referenceFastaFile |" : $referenceFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($chromosome = $1));
		$chromosomeSequenceLengthHash{$chromosome} += length($line);
	}
	close($reader);
}
my $pid2 = open2(my $reader2, my $writer2, "sort -t '\t' -k3,3 -k4,4 -k5,5n -k6,6n");
{
	my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2 -k3,3 -k4,4 -k5,5n -k6,6n");
	{
		open(my $reader, ($gtfFile =~ /\.gz$/) ? "gzip -dc $gtfFile |" : $gtfFile);
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^#/);
			my %tokenHash = ();
			@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line, -1);
			my %attributeHash = ();
			$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g);
			next unless($tokenHash{'feature'} eq 'exon');
			next unless(defined($chromosomeSequenceLengthHash{$tokenHash{'chromosome'}}));
			print $writer join("\t", @attributeHash{$geneAttribute, $transcriptAttribute}, @tokenHash{'chromosome', 'strand', 'start', 'end'}), "\n";
		}
		close($reader);
	}
	close($writer);
	{
		my ($geneTranscriptChromosomeStrand, $start, $end) = ('');
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{'gene', 'transcript', 'chromosome', 'strand', 'start', 'end'} = split(/\t/, $line, -1);
			if($geneTranscriptChromosomeStrand ne join("\t", @tokenHash{'gene', 'transcript', 'chromosome', 'strand'})) {
				print $writer2 join("\t", $geneTranscriptChromosomeStrand, $start, $end), "\n" if($geneTranscriptChromosomeStrand ne '');
				$geneTranscriptChromosomeStrand = join("\t", @tokenHash{'gene', 'transcript', 'chromosome', 'strand'});
				$start = $tokenHash{'start'};
			}
			$end = $tokenHash{'end'};
		}
		print $writer2 join("\t", $geneTranscriptChromosomeStrand, $start, $end), "\n" if($geneTranscriptChromosomeStrand ne '');
	}
	close($reader);
	waitpid($pid, 0);
}
close($writer2);
{
	my @forwardStartTokenHashList = ();
	my @forwardEndTokenHashList = ();
	my @reverseStartTokenHashList = ();
	my @reverseEndTokenHashList = ();
	while(my $line = <$reader2>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{'gene', 'transcript', 'chromosome', 'strand', 'start', 'end'} = split(/\t/, $line, -1);
		if($tokenHash{'strand'} eq '+') {
			printForward() if(@forwardStartTokenHashList && ($forwardStartTokenHashList[0]->{'chromosome'} ne $tokenHash{'chromosome'} || $forwardStartTokenHashList[0]->{'start'} < $tokenHash{'start'}));
			push(@forwardStartTokenHashList, \%tokenHash);
		}
		if($tokenHash{'strand'} eq '-') {
			printReverse() if(@reverseStartTokenHashList && ($reverseStartTokenHashList[0]->{'chromosome'} ne $tokenHash{'chromosome'} || $reverseStartTokenHashList[0]->{'start'} < $tokenHash{'start'}));
			push(@reverseStartTokenHashList, \%tokenHash);
		}
	}
	printForward() if(@forwardStartTokenHashList);
	printReverse() if(@reverseStartTokenHashList);
	printForward();
	printReverse();

	sub printForward {
		if(@forwardEndTokenHashList) {
			my ($chromosome, $start, $end) = ($forwardEndTokenHashList[0]->{'chromosome'}, $forwardEndTokenHashList[0]->{'end'} + 1, $forwardEndTokenHashList[0]->{'end'} + $maximumLength);
			$end = $chromosomeSequenceLengthHash{$chromosome} if($end > $chromosomeSequenceLengthHash{$chromosome});
			if($start <= $end) {
				unless(@forwardStartTokenHashList) {
					print join("\t", $chromosome, $start, $end, '+', $_->{'gene'}, $_->{'transcript'}), "\n" foreach(@forwardEndTokenHashList);
				} elsif($forwardEndTokenHashList[0]->{'chromosome'} ne $forwardStartTokenHashList[0]->{'chromosome'}) {
					print join("\t", $chromosome, $start, $end, '+', $_->{'gene'}, $_->{'transcript'}), "\n" foreach(@forwardEndTokenHashList);
				} elsif($forwardEndTokenHashList[0]->{'end'} < $forwardStartTokenHashList[0]->{'start'} - 1) {
					$end = $forwardStartTokenHashList[0]->{'start'} - 1 if($end > $forwardStartTokenHashList[0]->{'start'} - 1);
					print join("\t", $chromosome, $start, $end, '+', $_->{'gene'}, $_->{'transcript'}), "\n" foreach(@forwardEndTokenHashList);
				}
			}
		}
		foreach my $tokenHash (@forwardStartTokenHashList) {
			unless(@forwardEndTokenHashList) {
				push(@forwardEndTokenHashList, $tokenHash);
			} elsif($forwardEndTokenHashList[0]->{'chromosome'} ne $tokenHash->{'chromosome'}) {
				@forwardEndTokenHashList = ();
				push(@forwardEndTokenHashList, $tokenHash);
			} elsif($forwardEndTokenHashList[0]->{'end'} < $tokenHash->{'end'}) {
				@forwardEndTokenHashList = ();
				push(@forwardEndTokenHashList, $tokenHash);
			} elsif($forwardEndTokenHashList[0]->{'end'} == $tokenHash->{'end'}) {
				push(@forwardEndTokenHashList, $tokenHash);
			}
		}
		@forwardStartTokenHashList = ();
	}

	sub printReverse {
		if(@reverseStartTokenHashList) {
			my ($chromosome, $start, $end) = ($reverseStartTokenHashList[0]->{'chromosome'}, $reverseStartTokenHashList[0]->{'start'} - $maximumLength, $reverseStartTokenHashList[0]->{'start'} - 1);
			$start = 1 if($start < 1);
			if($start <= $end) {
				unless(@reverseEndTokenHashList) {
					print join("\t", $chromosome, $start, $end, '-', $_->{'gene'}, $_->{'transcript'}), "\n" foreach(@reverseStartTokenHashList);
				} elsif($reverseEndTokenHashList[0]->{'chromosome'} ne $reverseStartTokenHashList[0]->{'chromosome'}) {
					print join("\t", $chromosome, $start, $end, '-', $_->{'gene'}, $_->{'transcript'}), "\n" foreach(@reverseStartTokenHashList);
				} elsif($reverseEndTokenHashList[0]->{'end'} + 1 < $reverseStartTokenHashList[0]->{'start'}) {
					$start = $reverseEndTokenHashList[0]->{'end'} + 1 if($start < $reverseEndTokenHashList[0]->{'end'} + 1);
					print join("\t", $chromosome, $start, $end, '-', $_->{'gene'}, $_->{'transcript'}), "\n" foreach(@reverseStartTokenHashList);
				}
			}
		}
		foreach my $tokenHash (@reverseStartTokenHashList) {
			unless(@reverseEndTokenHashList) {
				push(@reverseEndTokenHashList, $tokenHash);
			} elsif($reverseEndTokenHashList[0]->{'chromosome'} ne $tokenHash->{'chromosome'}) {
				@reverseEndTokenHashList = ();
				push(@reverseEndTokenHashList, $tokenHash);
			} elsif($reverseEndTokenHashList[0]->{'end'} < $tokenHash->{'end'}) {
				@reverseEndTokenHashList = ();
				push(@reverseEndTokenHashList, $tokenHash);
			} elsif($reverseEndTokenHashList[0]->{'end'} == $tokenHash->{'end'}) {
				push(@reverseEndTokenHashList, $tokenHash);
			}
		}
		@reverseStartTokenHashList = ();
	}
}
close($reader2);
waitpid($pid2, 0);
