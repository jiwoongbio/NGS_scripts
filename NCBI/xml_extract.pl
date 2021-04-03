#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use XML::LibXML;

my ($xmlFile, $childNodeTagNames, @childNodeTagNamesList) = @ARGV;
my @childNodeTagNameList = split(/,/, $childNodeTagNames);
my @childNodeTagNameListList = map {[split(/,/, $_)]} @childNodeTagNamesList;
open(my $reader, $xmlFile);
binmode $reader;
my $dom = XML::LibXML->load_xml(IO => $reader);
my $root = $dom->documentElement();
if(@childNodeTagNameListList) {
	foreach my $childNode (getChildNodeList($root, @childNodeTagNameList)) {
		my @tokenList = ();
		foreach my $index (0 .. $#childNodeTagNameListList) {
			if(my @childNodeTagNameList = @{$childNodeTagNameListList[$index]}) {
				my $lastName = pop(@childNodeTagNameList);
				foreach my $childNode (getChildNodeList($childNode, @childNodeTagNameList)) {
					if(defined(my $attribute = $childNode->getAttribute($lastName))) {
						push(@{$tokenList[$index]}, $attribute);
					} else {
						push(@{$tokenList[$index]}, $_) foreach(map {$_->textContent ne '' ? $_->textContent : (map {$_->nodeName} $_->childNodes())} getChildNodeList($childNode, $lastName));
					}
				}
			} else {
				push(@{$tokenList[$index]}, $childNode->textContent);
			}
		}
		print join("\t", map {defined($_) ? join(',', @$_) : ''} @tokenList[0 .. $#childNodeTagNameListList]), "\n";
	}
} else {
	my $lastName = pop(@childNodeTagNameList);
	foreach my $childNode (getChildNodeList($root, @childNodeTagNameList)) {
		if(defined(my $attribute = $childNode->getAttribute($lastName))) {
			print $attribute, "\n";
		} else {
			print $_, "\n" foreach(map {$_->textContent ne '' ? $_->textContent : (map {$_->nodeName} $_->childNodes())} getChildNodeList($childNode, $lastName));
		}
	}
}
close($reader);

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
