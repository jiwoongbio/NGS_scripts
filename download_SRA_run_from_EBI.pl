#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'q' => \(my $quiet = ''),
	't=i' => \(my $maximumTryCount = 5),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   download_SRA_run_from_EBI.pl [options] SRA_run_accession [...]

Options: -h       display this help message
         -q       quiet
         -t INT   maximum try count [$maximumTryCount]

EOF
}
my $wgetOption = $quiet ? '--quiet' : '--no-verbose';
my (@runList) = @ARGV;
foreach my $run (@runList) {
	my @ftpList = ();
	my @md5List = ();
	open(my $reader, "wget $wgetOption -O - 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes' |");
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line, -1);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		push(@ftpList, $_) foreach(split(/;/, $tokenHash{'fastq_ftp'}));
		push(@md5List, $_) foreach(split(/;/, $tokenHash{'fastq_md5'}));
		die unless(scalar(@ftpList) == scalar(@md5List));
	}
	close($reader);
	foreach my $index (0 .. $#ftpList) {
		my ($ftp, $md5) = ($ftpList[$index], $md5List[$index]);
		(my $file = $ftp) =~ s/^.*\///;
		my $tryCount = 1;
		while($tryCount <= $maximumTryCount) {
			system("wget $wgetOption -O $file $ftp");
			chomp(my $md5sum = `md5sum $file`);
			$md5sum =~ s/\s+.*$//;
			last if($md5sum eq $md5);
			$tryCount += 1;
		}
		print STDERR "Failed: $ftp" if($tryCount > $maximumTryCount);
	}
}
