#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Time::HiRes;
use URI::Escape;
use XML::LibXML;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'e' => \(my $termEncoded = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   espell.pl [options] db term

Options: -h       display this help message
         -e       term encoded

EOF
}
my ($db, $term) = @ARGV;
my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $apiKey = $ENV{'NCBI_API_KEY'};
if($term eq '-') {
	while(my $term = <STDIN>) {
		chomp($term);
		print join("\t", $term, $_), "\n" foreach(espell($db, $term));
	}
} else {
	print join("\t", $term, $_), "\n" foreach(espell($db, $term));
}

sub espell {
	my ($db, $term) = @_;
	my $encodedTerm = $termEncoded ? $term : uri_escape($term);
	my $xmlString = '';
	while($xmlString !~ /<\/eSpellResult>\n?$/) {
		if(defined($apiKey)) {
			$xmlString = `wget --no-verbose -O - '$baseURL/espell.fcgi?db=$db&term=$encodedTerm&api_key=$apiKey'`;
			Time::HiRes::usleep(105);
		} else {
			$xmlString = `wget --no-verbose -O - '$baseURL/espell.fcgi?db=$db&term=$encodedTerm'`;
			Time::HiRes::usleep(350);
		}
	}
	my $dom = XML::LibXML->load_xml(string => $xmlString);
	my $root = $dom->documentElement();
	return map {$_->textContent} getChildNodeList($root, 'CorrectedQuery');
}

sub getChildNodeList {
	my ($node, @childNodeTagNameList) = @_;
	my @childNodeList = ();
	if(@childNodeTagNameList) {
		foreach my $childNode ($node->getChildrenByTagName(shift @childNodeTagNameList)) {
			push(@childNodeList, getChildNodeList($childNode, @childNodeTagNameList));
		}
	} else {
		push(@childNodeList, $node);
	}
	return @childNodeList;
}
