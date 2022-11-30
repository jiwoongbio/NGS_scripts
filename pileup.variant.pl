#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;
use List::Util qw(max);
use Getopt::Long qw(:config no_ignore_case);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'n=i' => \(my $numberPerThread = 1000),
	'd=i' => \(my $minimumReadBaseDepth = 0),
	'r=f' => \(my $minimumReadBaseDepthRatio = 0),
	'D=i' => \(my $minimumReadDepth = 0),
);
my $temporaryPrefix = "$temporaryDirectory/$hostname.$$";
system("rm -f $temporaryPrefix.*");
{
	my %pidHash = ();
	my @pidList = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			push(@pidList, $pid);
		} else {
			open($writer, "> $temporaryPrefix.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(@pidList) >= $number) {
			if($pidHash{$pidList[0]}) {
				my $pid = shift(@pidList);
				open(my $reader, "$temporaryPrefix.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryPrefix.$pid");
				delete $pidHash{$pid};
			}
			$pidHash{my $pid = wait()} = 1;
		}
	}
	sub forkPrint {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parentWriter)) {
			print $parentWriter @_;
		} else {
			print @_;
		}
	}
}
my ($pileupFile, $referenceFastaFile, @sampleList) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
open(my $reader, ($pileupFile =~ /\.gz$/ ? "gzip -dc $pileupFile |" : $pileupFile));
my @lineList = ();
while(my $line = <$reader>) {
	push(@lineList, $line);
	if(scalar(@lineList) >= $numberPerThread) {
		if($threads == 1) {
			printVariants(@lineList);
		} else {
			forkPrintSubroutine(\&printVariants, @lineList);
		}
		@lineList = ();
	}
}
if(@lineList) {
	if($threads == 1) {
		printVariants(@lineList);
	} else {
		forkPrintSubroutine(\&printVariants, @lineList);
	}
	@lineList = ();
}
forkPrintWait();
close($reader);

sub printVariants {
	foreach my $line (@_) {
		chomp($line);
		my @tokenList = split("\t", $line, -1);
		my ($chromosome, $position, $refBase) = @tokenList;
		$refBase =~ tr/a-z/A-Z/;
		my @readDepthList = ();
		my %readBaseDepthListHash = ();
		for(my $index = 3; $index < scalar(@tokenList); $index += 3) {
			my ($readDepth, $readBases, $baseQualities) = @tokenList[$index .. $index + 2];
			push(@readDepthList, $readDepth);
			foreach my $readBase (getReadBaseList($readBases)) {
				$readBase =~ s/^\^.//;
				$readBase =~ s/\$$//;
				$readBase =~ s/^[.,]/$refBase/;
				$readBase =~ tr/a-z/A-Z/;
				$readBaseDepthListHash{$readBase}->[($index / 3) - 1] += 1;
			}
		}
		next if(grep {$_ < $minimumReadDepth} @readDepthList);
		foreach my $readBase (grep {$_ ne $refBase && $_ ne '>' && $_ ne '<' && $_ !~ /^\*/} sort keys %readBaseDepthListHash) {
			my @readBaseDepthList = map {defined($_) ? $_ : 0} map {$readBaseDepthListHash{$readBase}->[$_]} 0 .. $#readDepthList;
			my @readBaseDepthRatioList = map {$readDepthList[$_] == 0 ? 0 : $readBaseDepthList[$_] / $readDepthList[$_]} 0 .. $#readDepthList;
			if(grep {$readBaseDepthList[$_] >= $minimumReadBaseDepth && $readBaseDepthRatioList[$_] >= $minimumReadBaseDepthRatio} 0 .. $#readDepthList) {
				printVariant($chromosome, $position, $refBase, $readBase, map {[$readDepthList[$_], $readBaseDepthList[$_], $readBaseDepthRatioList[$_]]} 0 .. $#readDepthList);
			}
		}
	}
}

sub printVariant {
	my ($chromosome, $position, $refBase, $altBase, @valueListList) = @_;
	while($altBase =~ s/([+-])([0-9]+)([A-Za-z]+)$//) {
		$refBase .= $3 if($1 eq '-');
		$altBase .= $3 if($1 eq '+');
	}
	($chromosome, $position, $refBase, $altBase) = extendIndel($chromosome, $position, $refBase, $altBase);
	if(@sampleList) {
		forkPrint(join("\t", $chromosome, $position, $refBase, $altBase, $sampleList[$_], @{$valueListList[$_]}), "\n") foreach(0 .. $#valueListList);
	} else {
		forkPrint(join("\t", $chromosome, $position, $refBase, $altBase, map {@$_} @valueListList), "\n");
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

sub extendIndel {
	my ($chromosome, $position, $refBase, $altBase) = @_;
	if($refBase ne $altBase) {
		while($refBase =~ /^$altBase/ || $altBase =~ /^$refBase/) {
			my $extBase = uc($db->seq($chromosome, ($_ = $position + length($refBase)), $_));
			$refBase = "$refBase$extBase";
			$altBase = "$altBase$extBase";
		}
		while(substr($refBase, -1, 1) eq substr($altBase, -1, 1) && !($refBase =~ /^$altBase/ || $altBase =~ /^$refBase/)) {
			substr($refBase, -1, 1, '');
			substr($altBase, -1, 1, '');
		}
		while($refBase =~ /$altBase$/ || $altBase =~ /$refBase$/) {
			my $extBase = uc($db->seq($chromosome, ($position = $position - 1), $position));
			$refBase = "$extBase$refBase";
			$altBase = "$extBase$altBase";
		}
		while(substr($refBase, 0, 1) eq substr($altBase, 0, 1) && !($refBase =~ /$altBase$/ || $altBase =~ /$refBase$/)) {
			$position = $position + 1;
			substr($refBase, 0, 1, '');
			substr($altBase, 0, 1, '');
		}
	}
	return ($chromosome, $position, $refBase, $altBase);
}
