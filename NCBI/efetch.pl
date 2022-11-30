#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

chomp(my $directory = `readlink -f $0`);
$directory =~ s/\/[^\/]*$//;

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage: efetch.pl [options] db id rettype retmode > record.txt

Options: -h       display this help message

Examples: efetch.pl nuccore NC_012920.1 gbwithparts text > NC_012920.1.gb
          efetch.pl pubmed 13054692 full xml > 13054692.xml

EOF
}

my ($db, $id, $rettype, $retmode) = @ARGV;
my $retmax = 100;
if($id eq '-') {
	chomp(my @idList = <STDIN>);
	for(my $index = 0; $index < scalar(@idList); $index += $retmax) {
		my $id = join(',', grep {defined} @idList[$index .. $index + $retmax - 1]);
		efetch($db, $id, $rettype, $retmode);
	}
} else {
	efetch($db, $id, $rettype, $retmode);
}

sub efetch {
	my ($db, $id, $rettype, $retmode) = @_;
	system("$directory/ecommon.pl 'efetch.fcgi?db=$db&id=$id&rettype=$rettype&retmode=$retmode'");
}
