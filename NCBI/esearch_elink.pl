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
	'mindate=s' => \(my $mindate = ''),
	'maxdate=s' => \(my $maxdate = ''),
	'l=s' => \(my $linkname = ''),
	't=s' => \(my $term = ''),
	'e' => \(my $termEncoded = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   esearch_elink.pl [options] dbfrom db termfrom

Options: -h       display this help message
         -mindate STR minimum date
         -maxdate STR maximum date
         -l STR   link name
         -t STR   term
         -e       term encoded

EOF
}
my ($dbfrom, $db, $termfrom) = @ARGV;
my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $apiKey = $ENV{'NCBI_API_KEY'};
my $retmax = 500;
if($term eq '-') {
	chomp(my @termList = <STDIN>);
	foreach my $term (@termList) {
		print join("\t", @$_), "\n" foreach(esearch_elink($dbfrom, $db, $termfrom, $mindate, $maxdate, $linkname, $term, $termEncoded));
	}
} else {
	print join("\t", @$_), "\n" foreach(esearch_elink($dbfrom, $db, $termfrom, $mindate, $maxdate, $linkname, $term, $termEncoded));
}

sub esearch_elink {
	my ($dbfrom, $db, $termfrom, $mindate, $maxdate, $linkname, $term, $termEncoded) = @_;
	my $encodedTermfrom = $termEncoded ? $termfrom : uri_escape($termfrom);
	$mindate = '' unless(defined($mindate));
	$maxdate = '' unless(defined($maxdate));
	$linkname = '' unless(defined($linkname));
	$term = '' unless(defined($term));
	my $encodedTerm = $termEncoded ? $term : uri_escape($term);
	my ($web, $key);
	{
		my $xmlString = '';
		while($xmlString !~ /<\/eSearchResult>\n?$/) {
			if(defined($apiKey)) {
				$xmlString = `wget --no-verbose -O - '$baseURL/esearch.fcgi?db=$dbfrom&term=$encodedTermfrom&usehistory=y&retmax=$retmax&mindate=$mindate&maxdate=$maxdate&api_key=$apiKey'`;
				Time::HiRes::usleep(105);
			} else {
				$xmlString = `wget --no-verbose -O - '$baseURL/esearch.fcgi?db=$dbfrom&term=$encodedTermfrom&usehistory=y&retmax=$retmax&mindate=$mindate&maxdate=$maxdate'`;
				Time::HiRes::usleep(350);
			}
		}
		my $dom = XML::LibXML->load_xml(string => $xmlString);
		my $root = $dom->documentElement();
		($key) = map {$_->textContent} getChildNodeList($root, 'QueryKey');
		($web) = map {$_->textContent} getChildNodeList($root, 'WebEnv');
	}
	my @idLinkNameList = ();
	{
		my $xmlString = '';
		while($xmlString !~ /<\/eLinkResult>\n?$/) {
			if(defined($apiKey)) {
				$xmlString = `wget --no-verbose -O - '$baseURL/elink.fcgi?dbfrom=$dbfrom&db=$db&query_key=$key&WebEnv=$web&linkname=$linkname&term=$encodedTerm&api_key=$apiKey'`;
				Time::HiRes::usleep(105);
			} else {
				$xmlString = `wget --no-verbose -O - '$baseURL/elink.fcgi?dbfrom=$dbfrom&db=$db&query_key=$key&WebEnv=$web&linkname=$linkname&term=$encodedTerm'`;
				Time::HiRes::usleep(350);
			}
		}
		my $dom = XML::LibXML->load_xml(string => $xmlString);
		my $root = $dom->documentElement();
		foreach my $linkSetDbNode (getChildNodeList($root, 'LinkSet', 'LinkSetDb')) {
			my @linkNameList = map {$_->textContent} getChildNodeList($linkSetDbNode, 'LinkName');
			push(@idLinkNameList, map {[$_->textContent, @linkNameList]} getChildNodeList($linkSetDbNode, 'Link', 'Id'));
		}
	}
	return @idLinkNameList;
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
