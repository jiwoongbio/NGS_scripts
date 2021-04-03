#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Time::HiRes;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   efetch.pl [options] db id rettype retmode > record.txt

Options: -h       display this help message

EOF
}
my ($db, $id, $rettype, $retmode) = @ARGV;
my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $apiKey = $ENV{'NCBI_API_KEY'};
my $retmax = 100;
if($id eq '-') {
	chomp(my @idList = <STDIN>);
	for(my $index = 0; $index < scalar(@idList); $index += $retmax) {
		my $id = join(',', grep {defined} @idList[$index .. $index + $retmax - 1]);
		print efetch($db, $id, $rettype, $retmode);
	}
} else {
	print efetch($db, $id, $rettype, $retmode);
}

sub efetch {
	my ($db, $id, $rettype, $retmode) = @_;
	my $output = '';
	while($output eq '') {
		if(defined($apiKey)) {
			$output = `wget --no-verbose -O - '$baseURL/efetch.fcgi?db=$db&id=$id&rettype=$rettype&retmode=$retmode&api_key=$apiKey'`;
			Time::HiRes::usleep(105);
		} else {
			$output = `wget --no-verbose -O - '$baseURL/efetch.fcgi?db=$db&id=$id&rettype=$rettype&retmode=$retmode'`;
			Time::HiRes::usleep(350);
		}
	}
	return $output;
}
