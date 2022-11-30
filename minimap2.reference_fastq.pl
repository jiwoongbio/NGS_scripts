#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
my %optionHash = ();
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'a' => \$optionHash{'a'},
	'Y' => \$optionHash{'Y'},
);
{
	my $parentPid = $$;
	my %pidHash = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			$pidHash{$pid} = 1;
		} else {
			open($writer, "> $temporaryDirectory/fork.$hostname.$parentPid.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(keys %pidHash) >= $number) {
			my $pid = wait();
			if($pidHash{$pid}) {
				open(my $reader, "$temporaryDirectory/fork.$hostname.$parentPid.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryDirectory/fork.$hostname.$parentPid.$pid");
				delete $pidHash{$pid};
			}
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
my @optionList = map {"-$_"} grep {$optionHash{$_}} keys %optionHash;
my ($referenceFastqFile, $fastqFile) = @ARGV;
{
	open(my $reader, ($referenceFastqFile =~ /\.gz$/ ? "gzip -dc $referenceFastqFile |" : $referenceFastqFile));
	while(my @lineList = getLineList($reader, 4)) {
		if($threads == 1) {
			printMapping(@lineList);
		} else {
			forkPrintSubroutine(\&printMapping, @lineList);
		}
	}
	forkPrintWait();
	close($reader);
}

sub printMapping {
	my $pid = open2(my $reader, my $writer, "minimap2 -x map-ont -t 1 @optionList - $fastqFile 2> /dev/null");
	print $writer "$_\n" foreach(@_);
	close($writer);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^@/);
		if($optionHash{'a'}) {
			my %tokenHash = ();
			(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
			$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
			next if($tokenHash{'flag'} & 4);
		}
		forkPrint("$line\n");
	}
	close($reader);
	waitpid($pid, 0);
}

sub getLineList {
	my ($reader, $number) = @_;
	my @lineList = ();
	foreach(1 .. $number) {
		if(defined(my $line = <$reader>)) {
			chomp($line);
			push(@lineList, $line);
		}
	}
	return @lineList;
}
