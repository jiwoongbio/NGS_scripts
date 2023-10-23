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
	'n=i' => \(my $numberInParallel = 0),
	's=i' => \(my $sleep = 1),
	'c' => \(my $doNotCanonicalize = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] script.sh [...]

Options: -h       display this help message
         -n INT   number of processes in parallel [$numberInParallel]
         -s INT   number of second to sleep [$sleep]
         -c       do not canonicalize

EOF
}
my (@shFileList) = @ARGV;
@shFileList = map {readlinkCanonicalize($_)} @shFileList if($doNotCanonicalize eq '');
chomp(my $nproc = `nproc`);
if($numberInParallel > 0) {
	exit if($numberInParallel > $nproc);
	my %pidReaderHash = ();
	while(1) {
		foreach my $pid (keys %pidReaderHash) {
			chomp(my $stat = `ps -p $pid --no-headers -o stat`);
			if($stat eq '' || $stat =~ /Z/) {
				close($pidReaderHash{$pid});
				delete($pidReaderHash{$pid});
			}
		}
		while(scalar(keys %pidReaderHash) < $numberInParallel && defined(my $shFile = shift(@shFileList))) {
			my $directory = $1 if($shFile =~ s/^(.*\/)//);
			my $pid = open(my $reader, "cd $directory; bash $shFile 1> $shFile.stdout 2> $shFile.stderr |");
			$pidReaderHash{$pid} = $reader;
		}
		last if(scalar(keys %pidReaderHash) == 0);
		sleep($sleep);
	}
} else {
	exit if(scalar(@shFileList) > $nproc);
	foreach my $shFile (@shFileList) {
		my $directory = $1 if($shFile =~ s/^(.*\/)//);
		system("cd $directory; bash $shFile 1> $shFile.stdout 2> $shFile.stderr &");
	}
}

sub readlinkCanonicalize {
	my @fileList = @_;
	chomp(@fileList = `readlink -f @fileList`);
	return @fileList;
}
