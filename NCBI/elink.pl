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
	'l=s' => \(my $linkname = ''),
	't=s' => \(my $term = ''),
	'e' => \(my $termEncoded = ''),
	'mindate=s' => \(my $mindate = ''),
	'maxdate=s' => \(my $maxdate = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   elink.pl [options] dbfrom db id

Options: -h       display this help message
         -l STR   link name
         -t STR   term
         -e       term encoded
         -mindate STR minimum date
         -maxdate STR maximum date

EOF
}
my ($dbfrom, $db, $id) = @ARGV;
my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $apiKey = $ENV{'NCBI_API_KEY'};
my $retmax = 100;
if($id eq '-') {
	chomp(my @idList = <STDIN>);
	for(my $index = 0; $index < scalar(@idList); $index += $retmax) {
		my $id = join(',', grep {defined} @idList[$index .. $index + $retmax - 1]);
		print join("\t", @$_), "\n" foreach(elink($dbfrom, $db, $id, $linkname, $term));
	}
} else {
	print join("\t", @$_), "\n" foreach(elink($dbfrom, $db, $id, $linkname, $term));
}

sub elink {
	my ($dbfrom, $db, $id, $linkname, $term) = @_;
	my $encodedTerm = $termEncoded ? $term : uri_escape($term);
	my @idLinkNameList = ();
	my $xmlString = '';
	while($xmlString !~ /<\/eLinkResult>\n?$/) {
		if(defined($apiKey)) {
			$xmlString = `wget --no-verbose -O - '$baseURL/elink.fcgi?dbfrom=$dbfrom&db=$db&id=$id&linkname=$linkname&term=$encodedTerm&mindate=$mindate&maxdate=$maxdate&api_key=$apiKey'`;
			Time::HiRes::usleep(105);
		} else {
			$xmlString = `wget --no-verbose -O - '$baseURL/elink.fcgi?dbfrom=$dbfrom&db=$db&id=$id&linkname=$linkname&term=$encodedTerm&mindate=$mindate&maxdate=$maxdate'`;
			Time::HiRes::usleep(350);
		}
	}
	my $dom = XML::LibXML->load_xml(string => $xmlString);
	my $root = $dom->documentElement();
	foreach my $linkSetDbNode (getChildNodeList($root, 'LinkSet', 'LinkSetDb')) {
		my @linkNameList = map {$_->textContent} getChildNodeList($linkSetDbNode, 'LinkName');
		push(@idLinkNameList, map {[$_->textContent, @linkNameList]} getChildNodeList($linkSetDbNode, 'Link', 'Id'));
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
