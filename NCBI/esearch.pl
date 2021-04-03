#!/bin/env perl
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
	'mindate=s' => \(my $mindate = ''),
	'maxdate=s' => \(my $maxdate = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   esearch.pl [options] db term > id.txt

Options: -h       display this help message
         -e       term encoded
         -mindate STR minimum date
         -maxdate STR maximum date

EOF
}
my ($db, $term) = @ARGV;
my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $apiKey = $ENV{'NCBI_API_KEY'};
my $retmax = 500;
if($term eq '-') {
	chomp(my @termList = <STDIN>);
	foreach my $term (@termList) {
		print "$_\n" foreach(esearch($db, $term));
	}
} else {
	print "$_\n" foreach(esearch($db, $term));
}

sub esearch {
	my ($db, $term) = @_;
	my $encodedTerm = $termEncoded ? $term : uri_escape($term);
	my @idList = ();
	for(my $retstart = 0, my ($count, $web, $key) = ($retmax); $retstart < $count;) {
		my $xmlString = '';
		unless(defined($web) && defined($key)) {
			if(defined($apiKey)) {
				$xmlString = `wget --no-verbose -O - '$baseURL/esearch.fcgi?db=$db&term=$encodedTerm&usehistory=y&retmax=$retmax&retstart=$retstart&mindate=$mindate&maxdate=$maxdate&api_key=$apiKey'`;
			} else {
				$xmlString = `wget --no-verbose -O - '$baseURL/esearch.fcgi?db=$db&term=$encodedTerm&usehistory=y&retmax=$retmax&retstart=$retstart&mindate=$mindate&maxdate=$maxdate'`;
			}
		} else {
			if(defined($apiKey)) {
				$xmlString = `wget --no-verbose -O - '$baseURL/esearch.fcgi?WebEnv=$web&query_key=$key&retmax=$retmax&retstart=$retstart&api_key=$apiKey'`;
			} else {
				$xmlString = `wget --no-verbose -O - '$baseURL/esearch.fcgi?WebEnv=$web&query_key=$key&retmax=$retmax&retstart=$retstart'`;
			}
		}
		if($xmlString =~ /<\/eSearchResult>\n?$/) {
			my $dom = XML::LibXML->load_xml(string => $xmlString);
			my $root = $dom->documentElement();
			($count) = map {$_->textContent} getChildNodeList($root, 'Count');
			unless(defined($web) && defined($key)) {
				($web) = map {$_->textContent} getChildNodeList($root, 'WebEnv');
				($key) = map {$_->textContent} getChildNodeList($root, 'QueryKey');
			}
			push(@idList, map {$_->textContent} getChildNodeList($root, 'IdList', 'Id'));
			$retstart += $retmax;
		}
		if(defined($apiKey)) {
			Time::HiRes::usleep(105);
		} else {
			Time::HiRes::usleep(350);
		}
	}
	return @idList;
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
