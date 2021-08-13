#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
my @coverageDepthList = ();
GetOptions(
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 1000),
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'S=s' => \(my $stranded = ''),
	'T' => \(my $printTotal = ''),
	'c=i' => \@coverageDepthList,
	'doNotUseSamModule' => \(my $doNotUseSamModule = ''),
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
my $temporaryPrefix = "$temporaryDirectory/MetaPrism.$hostname.$$";
my ($regionFile, @bamFileList) = @ARGV;
my $bamFile = '';
my $samModule = $doNotUseSamModule ? '' : eval {require Bio::DB::Sam; 1;};
if($samModule) {
	foreach my $bamFile (@bamFileList) {
		unless(-r "$bamFile.bai") {
			(my $baiFile = $bamFile) =~ s/\.bam/.bai/;
			unless(-r $baiFile) {
				system("samtools index $bamFile");
			} else {
				system("ln --relative --symbolic $baiFile $bamFile.bai");
			}
		}
	}
} else {
	if(scalar(@bamFileList) == 1) {
		($bamFile) = @bamFileList;
		unless(-r "$bamFile.bai") {
			(my $baiFile = $bamFile) =~ s/\.bam/.bai/;
			unless(-r $baiFile) {
				system("samtools index $bamFile");
			}
		}
	} elsif(scalar(@bamFileList) > 1) {
		$bamFile = "$temporaryPrefix.bam";
		system("samtools merge --threads $threads $bamFile @bamFileList");
		system("samtools index $bamFile");
	}
}
my $pid = open2(my $reader, my $writer, 'sort -n');
forkPrintParentWriter($writer);
{
	open(my $reader, $regionFile);
	my $order = 0;
	my @tokenListList = ();
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		push(@tokenListList, [$order += 1, @tokenList]);
		if(scalar(@tokenListList) >= $numberPerThread) {
			if($threads == 1) {
				printDepth(@tokenListList);
			} else {
				forkPrintSubroutine(\&printDepth, @tokenListList);
			}
			@tokenListList = ();
		}
	}
	if(@tokenListList) {
		if($threads == 1) {
			printDepth(@tokenListList);
		} else {
			forkPrintSubroutine(\&printDepth, @tokenListList);
		}
	}
	close($reader);
	forkPrintWait();
}
forkPrintParentWriter();
close($writer);
{
	my $totalLength = 0;
	my $totalDepth = 0;
	my @totalCoverageList = ();
	while(my $line = <$reader>) {
		chomp($line);
		my ($order, @tokenList) = split(/\t/, $line, -1);
		my ($length, $depth, @coverageList) = split(/,/, pop(@tokenList));
		print join("\t", @tokenList, $depth / $length, map {sprintf('%.1f%%', $_ / $length * 100)} @coverageList), "\n" if($printTotal eq '');
		$totalLength += $length;
		$totalDepth += $depth;
		$totalCoverageList[$_] += $coverageList[$_] foreach(0 .. $#coverageList);
	}
	if($printTotal) {
		print join("\t", 'Mean depth', $totalDepth / $totalLength), "\n";
		print join("\t", "% Depth >= $coverageDepthList[$_]", sprintf('%.1f%%', $totalCoverageList[$_] / $totalLength * 100)), "\n" foreach(0 .. $#totalCoverageList);
	}
}
close($reader);
waitpid($pid, 0);
unless($samModule) {
	if(scalar(@bamFileList) > 1) {
		$bamFile = "$temporaryPrefix.bam";
		system("rm $bamFile $bamFile.bai");
	}
}

sub printDepth {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} @bamFileList if($samModule);
	foreach(@_) {
		my ($order, @tokenList) = @$_;
		my ($chromosome, $start, $end, $strand) = @tokenList;
		my $depth = 0;
		my @coverageList = (0) x scalar(@coverageDepthList);
		if($samModule) {
			my %positionDepthHash = ();
			foreach my $sam (@samList) {
				foreach my $alignment ($sam->get_features_by_location(-seq_id => $chromosome, -start => $start, -end => $end)) {
					next if($alignment->qual < $minimumMappingQuality);
					next if(($alignment->flag & int($includeFlag)) != $includeFlag);
					next if($alignment->flag & int($excludeFlag));
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
						next unless(grep {($alignment->flag & 253) == $_} 97, 145, 0);
					}
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
						next unless(grep {($alignment->flag & 253) == $_} 81, 161, 16);
					}
					my @positionList = getPositionList($alignment->start, $alignment->cigar_str);
					@positionList = grep {$_ ne ''} @positionList;
					@positionList = grep {$start <= $_ && $_ <= $end} @positionList;
					$positionDepthHash{$_} += 1 foreach(@positionList);
				}
			}
			foreach my $positionDepth (values %positionDepthHash) {
				$depth += $positionDepth;
				$coverageList[$_] += 1 foreach(grep {$positionDepth >= $coverageDepthList[$_]} 0 .. $#coverageDepthList);
			}
		} else {
			if($bamFile ne '') {
				my $reader;
				if($stranded eq '') {
					open($reader, "samtools view -h -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$start-$end | samtools depth -d 0 - |");
				} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
					open($reader, "samtools view -h -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$start-$end | perl -Mstrict -Mwarnings -wne 'chomp(my \$line = \$_); my \@tokenList = split(/\\t/, \$line, -1); print \"\$line\\n\" if(\$line =~ /^\@/ || grep {(\$tokenList[1] & 253) == \$_} 97, 145, 0);' | samtools depth -d 0 - |");
				} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
					open($reader, "samtools view -h -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$start-$end | perl -Mstrict -Mwarnings -wne 'chomp(my \$line = \$_); my \@tokenList = split(/\\t/, \$line, -1); print \"\$line\\n\" if(\$line =~ /^\@/ || grep {(\$tokenList[1] & 253) == \$_} 81, 161, 16);' | samtools depth -d 0 - |");
				}
				while(my $line = <$reader>) {
					chomp($line);
					my ($chromosome, $position, $positionDepth) = split(/\t/, $line, -1);
					if($start <= $position && $position <= $end) {
						$depth += $positionDepth;
						$coverageList[$_] += 1 foreach(grep {$positionDepth >= $coverageDepthList[$_]} 0 .. $#coverageDepthList);
					}
				}
				close($reader);
			}
		}
		my $length = $end - $start + 1;
		forkPrint(join("\t", $order, @tokenList, join(',', $length, $depth, @coverageList)), "\n");
	}
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		} elsif($operation eq 'D') {
			$position += $length;
		} elsif($operation eq 'N') {
			$position += $length;
		} elsif($operation eq 'S') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		}
	}
	return @positionList;
}
