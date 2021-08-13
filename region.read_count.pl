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
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 1000),
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'S=s' => \(my $stranded = ''),
	'T' => \(my $printTotal = ''),
	'c' => \(my $countFragment = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl region.read_count.pl [options] region.txt sample.sorted.bam [...] > region.read_count.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -f INT   include flag [$includeFlag]
         -F INT   exclude flag [$excludeFlag]
         -S STR   stranded, "f" or "r"
         -T       print total read count
         -c       count fragment

EOF
}
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
my @samMandatoryColumnList = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');
my $samModule = eval {require Bio::DB::Sam; 1;};
my ($regionFile, @bamFileList) = @ARGV;
foreach my $bamFile (@bamFileList) {
	if(-r "$bamFile.bai") {
	} else {
		(my $baiFile = $bamFile) =~ s/\.bam/.bai/;
		if(-r $baiFile) {
			system("ln -s $baiFile $bamFile.bai");
		} else {
			system("samtools index $bamFile");
		}
	}
}
my $pid = open2(my $reader, my $writer, $printTotal ? 'sort -u | wc -l' : "sort -t '\t' -k1,1n -k2,2 | uniq | cut -f3- | uniq -c");
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
				printRead(@tokenListList);
			} else {
				forkPrintSubroutine(\&printRead, @tokenListList);
			}
			@tokenListList = ();
		}
	}
	if(@tokenListList) {
		if($threads == 1) {
			printRead(@tokenListList);
		} else {
			forkPrintSubroutine(\&printRead, @tokenListList);
		}
	}
	close($reader);
	forkPrintWait();
}
forkPrintParentWriter();
close($writer);
if($printTotal) {
	while(my $line = <$reader>) {
		chomp($line);
		print "$line\n";
	}
} else {
	while(my $line = <$reader>) {
		chomp($line);
		$line = join("\t", $line, $1) if($line =~ s/^ *([0-9]+) //);
		print "$line\n";
	}
}
close($reader);
waitpid($pid, 0);

sub printRead {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} @bamFileList if($samModule);
	foreach(@_) {
		my ($order, @tokenList) = @$_;
		my ($chromosome, $start, $end, $strand) = @tokenList;
		if($samModule) {
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
					my $readName = $countFragment ? $alignment->query->name : join('/', $alignment->query->name, ($alignment->flag & 192) / 64);
					if($printTotal) {
						forkPrint(join("\t", $readName), "\n");
					} else {
						forkPrint(join("\t", $order, $readName, @tokenList), "\n");
					}
				}
			}
		} else {
			foreach my $bamFile (@bamFileList) {
				open(my $reader, "samtools view -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$start-$end |");
				while(my $line = <$reader>) {
					chomp($line);
					my %tokenHash = ();
					(@tokenHash{@samMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line);
					$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
						next unless(grep {($tokenHash{'flag'} & 253) == $_} 97, 145, 0);
					}
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
						next unless(grep {($tokenHash{'flag'} & 253) == $_} 81, 161, 16);
					}
					my $readName = $countFragment ? $tokenHash{'qname'} : join('/', $tokenHash{'qname'}, ($tokenHash{'flag'} & 192) / 64);
					if($printTotal) {
						forkPrint(join("\t", $readName), "\n");
					} else {
						forkPrint(join("\t", $order, $readName, @tokenList), "\n");
					}
				}
				close($reader);
			}
		}
	}
}
