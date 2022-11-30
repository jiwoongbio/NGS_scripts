#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'd=s' => \(my $datafile = 'EDNAFULL'),
	'O=f' => \(my $gapOpen = ''),
	'E=f' => \(my $gapExtend = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   alignment.position_base.pl [options] reference.fasta alternate.fasta > position_base.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -d STR   scoring matrix
         -O FLOAT gap open penalty
         -E FLOAT gap extend penalty

EOF
}
my ($referenceFastaFile, $alternateFastaFile) = @ARGV;
my @referenceNameList = ();
my %referenceNameSequenceHash = ();
{
	open(my $reader, ($referenceFastaFile =~ /\.gz$/ ? "gzip -dc $referenceFastaFile |" : $referenceFastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		next if($line eq '');
		if($line =~ /^>(\S*)/) {
			push(@referenceNameList, $1);
		} else {
			$referenceNameSequenceHash{$referenceNameList[-1]} .= uc($line);
		}
	}
	close($reader);
}
my @alternateNameList = ();
my %alternateNameSequenceHash = ();
{
	open(my $reader, ($alternateFastaFile =~ /\.gz$/ ? "gzip -dc $alternateFastaFile |" : $alternateFastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		next if($line eq '');
		if($line =~ /^>(\S*)/) {
			push(@alternateNameList, $1);
		} else {
			$alternateNameSequenceHash{$alternateNameList[-1]} .= uc($line);
		}
	}
	close($reader);
}
my @alignmentIdentityList = ();
foreach my $referenceName (@referenceNameList) {
	foreach my $alternateName (@alternateNameList) {
		my $referenceSequence = $referenceNameSequenceHash{$referenceName};
		my $alternateSequence = $alternateNameSequenceHash{$alternateName};
		my ($referenceAlignment, $alternateAlignment, $identity) = ($referenceSequence, $alternateSequence, 1);
		if($referenceAlignment ne $alternateAlignment) {
			if($alternateAlignment =~ /$referenceAlignment/) {
				$identity = length($referenceAlignment) / length($alternateAlignment);
				$referenceAlignment = ('-' x index($alternateAlignment, $referenceAlignment)) . $referenceAlignment;
				$referenceAlignment .= '-' x length($1) if(substr($alternateAlignment, length($referenceAlignment)) =~ /^(.*[^A])A*$/);
			} elsif($referenceAlignment =~ /$alternateAlignment/) {
				$identity = length($alternateAlignment) / length($referenceAlignment);
				$alternateAlignment = ('-' x index($referenceAlignment, $alternateAlignment)) . $alternateAlignment;
				$alternateAlignment .= '-' x (length($referenceAlignment) - length($alternateAlignment));
			} else {
				($referenceAlignment, $alternateAlignment, $identity) = needle($referenceSequence, $alternateSequence, $datafile, $gapOpen, $gapExtend);
				($referenceAlignment, $alternateAlignment, $identity) = stretcher($referenceSequence, $alternateSequence, $datafile, $gapOpen, $gapExtend) if($identity eq '');
			}
		}
		push(@alignmentIdentityList, [$referenceName, $referenceAlignment, $alternateName, $alternateAlignment, $identity]);
	}
}
@alignmentIdentityList = sort {$b->[-1] <=> $a->[-1]} @alignmentIdentityList;
my %referenceNameHash = ();
my %alternateNameHash = ();
foreach(@alignmentIdentityList) {
	my ($referenceName, $referenceAlignment, $alternateName, $alternateAlignment, $identity) = @$_;
	my @referenceBaseList = split(//, $referenceAlignment);
	my @alternateBaseList = split(//, $alternateAlignment);
	my $referencePosition = 0;
	my $alternatePosition = 0;
	foreach my $index (0 .. $#referenceBaseList) {
		my $referenceBase = $referenceBaseList[$index];
		my $alternateBase = $alternateBaseList[$index];
		$referencePosition += 1 if($referenceBase ne '-');
		$alternatePosition += 1 if($alternateBase ne '-');
		print join("\t", $referenceName, $referencePosition, $referenceBase, $alternateName, $alternatePosition, $alternateBase), "\n";
	}
	$referenceNameHash{$referenceName} = 1;
	$alternateNameHash{$alternateName} = 1;
	last if(scalar(keys %referenceNameHash) == scalar(@referenceNameList));
	last if(scalar(keys %alternateNameHash) == scalar(@alternateNameList));
}

sub needle {
	my ($seqA, $seqB, $datafile, $gapOpen, $gapExtend) = @_;
	$gapOpen = '' unless(defined($gapOpen));
	$gapExtend = '' unless(defined($gapExtend));
	my $prefix = "$temporaryDirectory/needle.$hostname.$$";
	system("rm -fr $prefix.*");
	{
		open(my $writer, "> $prefix.asequence");
		print $writer ">seqA\n";
		print $writer "$seqA\n";
		close($writer);
	}
	{
		open(my $writer, "> $prefix.bsequence");
		print $writer ">seqB\n";
		print $writer "$seqB\n";
		close($writer);
	}
	print STDERR `echo -en '$gapOpen\\n$gapExtend\\n' | timeout 600 needle -asequence $prefix.asequence -bsequence $prefix.bsequence -outfile $prefix.outfile -datafile $datafile`;
	print STDERR "\n";
	my ($alignA, $alignB, $identity) = ('', '', '');
	if(-e "$prefix.outfile") {
		open(my $reader, "$prefix.outfile");
		while(my $line = <$reader>) {
			chomp($line);
			$identity = $1 / $2 if($line =~ /^# Identity: +([0-9]+)\/([0-9]+) \( *[0-9.]+%\)$/);
			$alignA .= $1 if($line =~ /^seqA +[0-9]+ +([^ ]+) +[0-9]+$/);
			$alignB .= $1 if($line =~ /^seqB +[0-9]+ +([^ ]+) +[0-9]+$/);
		}
		close($reader);
	}
	system("rm -fr $prefix.*");
	return ($alignA, $alignB, $identity);
}

sub stretcher {
	my ($seqA, $seqB, $datafile, $gapOpen, $gapExtend) = @_;
	$gapOpen = defined($gapOpen) && $gapOpen ne '' ? "-gapopen $gapOpen" : '';
	$gapExtend = defined($gapExtend) && $gapExtend ne '' ? "-gapextend $gapExtend" : '';
	my $prefix = "$temporaryDirectory/stretcher.$hostname.$$";
	system("rm -fr $prefix.*");
	{
		open(my $writer, "> $prefix.asequence");
		print $writer ">seqA\n";
		print $writer "$seqA\n";
		close($writer);
	}
	{
		open(my $writer, "> $prefix.bsequence");
		print $writer ">seqB\n";
		print $writer "$seqB\n";
		close($writer);
	}
	print STDERR `timeout 600 stretcher -asequence $prefix.asequence -bsequence $prefix.bsequence -outfile $prefix.outfile -datafile $datafile $gapOpen $gapExtend`;
	print STDERR "\n";
	my ($alignA, $alignB, $identity) = ('', '', '');
	if(-e "$prefix.outfile") {
		open(my $reader, "$prefix.outfile");
		while(my $line = <$reader>) {
			chomp($line);
			$identity = $1 / $2 if($line =~ /^# Identity: +([0-9]+)\/([0-9]+) \( *[0-9.]+%\)$/);
			$alignA .= $1 if($line =~ /^ *seqA (.+)$/);
			$alignB .= $1 if($line =~ /^ *seqB (.+)$/);
		}
		close($reader);
	}
	system("rm -fr $prefix.*");
	return ($alignA, $alignB, $identity);
}
