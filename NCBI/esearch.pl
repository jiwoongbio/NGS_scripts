#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

chomp(my $directory = `readlink -f $0`);
$directory =~ s/\/[^\/]*$//;

use URI::Escape;
use XML::LibXML;
use Getopt::Long qw(:config no_ignore_case);

my @defaultQueryStringList = ();
GetOptions(
	'h' => \(my $help = ''),
	'e' => \(my $termEncoded = ''),
	'mindate=s' => \(my $mindate = ''),
	'maxdate=s' => \(my $maxdate = ''),
	'q=s' => \@defaultQueryStringList,
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage: esearch.pl [options] db term > id.txt

Options: -h       display this help message
         -e       term encoded
         -mindate STR minimum date
         -maxdate STR maximum date
         -q       default query string

Examples: esearch.pl nuccore 'txid9606[Organism] AND biomol_rRNA[Properties]' > human.rRNA.nuccore.txt
          esearch.pl biosample 'antibiogram[filter]' > antibiogram.biosample.txt

EOF
}

my ($db, $term) = @ARGV;
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
		my @queryStringList = ();
		unless(defined($web) && defined($key)) {
			push(@queryStringList, "db=$db&term=$encodedTerm&usehistory=y&mindate=$mindate&maxdate=$maxdate");
		} else {
			push(@queryStringList, "WebEnv=$web&query_key=$key");
		}
		push(@queryStringList, "retmax=$retmax&retstart=$retstart");
		my $queryString = join('&', @queryStringList, @defaultQueryStringList);
		my $xmlString = `$directory/ecommon.pl 'esearch.fcgi?$queryString'`;
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
