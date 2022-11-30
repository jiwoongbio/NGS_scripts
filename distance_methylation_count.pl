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
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 1000),
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
my ($tssFile, $bgzMethylationFile, $maximumDistance) = @ARGV;
{
	open(my $reader, $tssFile);
	my @lineList = ();
	while(my $line = <$reader>) {
		push(@lineList, $line);
		if(scalar(@lineList) >= $numberPerThread) {
			if($threads == 1) {
				printDistanceMethylationCount(@lineList);
			} else {
				forkPrintSubroutine(\&printDistanceMethylationCount, @lineList);
			}
			@lineList = ();
		}
	}
	if(@lineList) {
		if($threads == 1) {
			printDistanceMethylationCount(@lineList);
		} else {
			forkPrintSubroutine(\&printDistanceMethylationCount, @lineList);
		}
	}
	close($reader);
	forkPrintWait();
}

sub printDistanceMethylationCount {
	chomp(my @lineList = @_);
	foreach my $line (@lineList) {
		my @tokenList = split(/\t/, $line, -1);
		my ($chromosome, $position, $strand) = @tokenList;
		my ($start, $end) = ('', '');
		if($strand eq '+') {
			$start = $position - $maximumDistance;
			$end = $position + ($maximumDistance - 1);
		}
		if($strand eq '-') {
			$start = $position - ($maximumDistance - 1);
			$end = $position + $maximumDistance;
		}
		my %distanceMethylationCountHash = ();
		open(my $reader, "tabix $bgzMethylationFile $chromosome:$start-$end |");
		while(my $line = <$reader>) {
			chomp($line);
			my ($chromosome, $distance, $methylation) = split(/\t/, $line, -1);
			$distance = $distance - $position if($strand eq '+');
			$distance = $position - $distance if($strand eq '-');
			$distance += 1 if($distance >= 0);
			$distanceMethylationCountHash{$distance}->{$methylation} += 1;
		}
		close($reader);
		foreach my $distance (sort {$a <=> $b} keys %distanceMethylationCountHash) {
			foreach my $methylation (sort keys %{$distanceMethylationCountHash{$distance}}) {
				my $count = $distanceMethylationCountHash{$distance}->{$methylation};
				forkPrint(join("\t", @tokenList, $distance, $methylation, $count), "\n");
			}
		}
	}
}
