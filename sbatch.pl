#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
	'p=s' => \(my $partition = '128GB'),
	'g=i' => \(my $memoryGB = 0),
	'n=i' => \(my $numberInParallel = 1),
	'm=i' => \(my $numberInNode = 1),
	'c' => \(my $doNotCanonicalize = ''),
	'f=s' => \(my $shFileFile = ''),
	'w=i' => \(my $wait = 0),
	'a' => \(my $addToRunningNode = ''),
);
my @partitionMemoryGBList = (['32GB', 30], ['128GB', 124], ['256GB', 250]);
my $partitionMemoryGBs = join(', ', map {"$_->[1] for $_->[0]"} @partitionMemoryGBList);
chomp(my $scriptDirectory = `readlink -f $0`);
$scriptDirectory =~ s/\/[^\/]*$//;
my @runningNodeList = ();
if($addToRunningNode) {
	open(my $reader, "perl $scriptDirectory/squeue.pl -u $ENV{'USER'} |");
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		push(@runningNodeList, $tokenList[7]) if($tokenList[1] eq $partition && $tokenList[4] eq 'RUNNING');
	}
	close($reader);
}
$memoryGB = $_->[1] foreach(grep {$_->[0] eq $partition} @partitionMemoryGBList);
my $memoryMB = $memoryGB * 1024;
my (@shFileList) = @ARGV;
if($shFileFile ne '') {
	open(my $reader, $shFileFile);
	while(my $line = <$reader>) {
		chomp($line);
		push(@shFileList, split(/\s+/, $line));
	}
	close($reader);
}
if($help || scalar(@shFileList) == 0) {
	die <<EOF;

Usage:   $file [options] script.sh [...]

Options: -h       display this help message
         -p STR   partition [$partition]
         -g INT   memory in GB [$partitionMemoryGBs]
         -n INT   number of processes in parallel [$numberInParallel]
         -m INT   number of processes in node [$numberInNode]
         -c       do not canonicalize

EOF
}
@shFileList = map {readlinkCanonicalize($_)} @shFileList if($doNotCanonicalize eq '');
if($numberInParallel == 1 && $numberInNode == 1) {
	foreach my $shFile (@shFileList) {
		sleep(10) while($wait > 0 && scalar(grep {$_->{'UserId'} =~ /^$ENV{'USER'}\([0-9]*\)$/ && $_->{'Partition'} eq $partition} getTokenHashList()) >= $wait);
		if($addToRunningNode) {
			if(@runningNodeList) {
				my $runningNode = shift(@runningNodeList);
				system("ssh $runningNode 'perl $scriptDirectory/bash.pl -c $shFile'");
			} else {
				print STDERR "$shFile\n";
			}
			next;
		}
		my $directory = $1 if($shFile =~ s/^(.*\/)//);
		system("cd $directory; sbatch --partition=$partition --mem $memoryMB --time=14-00:00:00 --workdir=\$PWD --job-name=$shFile --output=$shFile.stdout --error=$shFile.stderr $shFile");
		sleep(10);
	}
} elsif($numberInParallel >= $numberInNode) {
	@shFileList = grep {$_ !~ /\/sbatch[0-9]{3}\.sh$/} @shFileList;
	my $number = 1;
	while(my @shFileList = splice(@shFileList, 0, $numberInParallel)) {
		if($addToRunningNode) {
			if(@runningNodeList) {
				my $runningNode = shift(@runningNodeList);
				system("ssh $runningNode 'perl $scriptDirectory/bash.pl -c @shFileList'");
			} else {
				print STDERR "$_\n" foreach(@shFileList);
			}
			next;
		}
		my $shFile;
		$number += 1 while(-e ($shFile = sprintf('sbatch%03d.sh', $number)));
		open(my $writer, "> $shFile");
		print $writer "#!/bin/bash\n";
		print $writer "rm -rf /tmp/* 2> /dev/null\n";
		print $writer "rm -rf /dev/shm/* 2> /dev/null\n";
		foreach my $shFile (@shFileList) {
			my $directory = $1 if($shFile =~ s/^(.*\/)//);
			print $writer "cd $directory; sh $shFile 1> $shFile.stdout 2> $shFile.stderr &\n";
		}
		print $writer "wait\n";
		close($writer);
		sleep(10) while($wait > 0 && scalar(grep {$_->{'UserId'} =~ /^$ENV{'USER'}\([0-9]*\)$/ && $_->{'Partition'} eq $partition} getTokenHashList()) >= $wait);
		system("sbatch --partition=$partition --mem $memoryMB --time=14-00:00:00 --job-name=$shFile --output=$shFile.stdout --error=$shFile.stderr $shFile");
		sleep(10);
	}
} elsif($numberInParallel < $numberInNode) {
	@shFileList = grep {$_ !~ /\/sbatch[0-9]{3}\.sh$/} @shFileList;
	my $number = 1;
	while(my @shFileList = splice(@shFileList, 0, $numberInNode)) {
		if($addToRunningNode) {
			if(@runningNodeList) {
				my $runningNode = shift(@runningNodeList);
				system("ssh $runningNode 'perl $scriptDirectory/bash.pl -n $numberInParallel -c @shFileList'");
			} else {
				print STDERR "$_\n" foreach(@shFileList);
			}
			next;
		}
		my $shFile;
		$number += 1 while(-e ($shFile = sprintf('sbatch%03d.sh', $number)));
		open(my $writer, "> $shFile");
		print $writer "#!/bin/bash\n";
		print $writer "rm -rf /tmp/* 2> /dev/null\n";
		print $writer "rm -rf /dev/shm/* 2> /dev/null\n";
		print $writer "perl $scriptDirectory/bash.pl -n $numberInParallel -c @shFileList\n";
		close($writer);
		sleep(10) while($wait > 0 && scalar(grep {$_->{'UserId'} =~ /^$ENV{'USER'}\([0-9]*\)$/ && $_->{'Partition'} eq $partition} getTokenHashList()) >= $wait);
		system("sbatch --partition=$partition --mem $memoryMB --time=14-00:00:00 --job-name=$shFile --output=$shFile.stdout --error=$shFile.stderr $shFile");
		sleep(10);
	}
}

sub getTokenHashList {
	open(my $reader, 'scontrol show job |');
	my @tokenHashList = ();
	my $tokenHash = {};
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^JobId=/) {
			push(@tokenHashList, $tokenHash) if(%$tokenHash);
			$tokenHash = {};
		}
		$tokenHash->{$_->[0]} = $_->[1] foreach(grep {scalar(@$_) == 2} map {[split(/=/, $_, 2)]} split(/\s+/, $line, -1));
	}
	push(@tokenHashList, $tokenHash) if(%$tokenHash);
	close($reader);
	return @tokenHashList;
}

sub readlinkCanonicalize {
	my @fileList = @_;
	chomp(@fileList = `readlink -f @fileList`);
	return @fileList;
}
