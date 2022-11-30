#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my @unwrapTagList = ();
GetOptions(
	'u=s' => \@unwrapTagList,
);
@unwrapTagList = map {split(/,/, $_)} @unwrapTagList;
my ($htmlFile) = @ARGV;
my $html = '';
{
	open(my $reader, $htmlFile);
	$html .= $_ while(<$reader>);
	close($reader);
}
$html = unwrapTag($html, $_) foreach(@unwrapTagList);
$html =~ s/&nbsp;/ /g;
$html =~ s/&#x000a0;/ /g;
$html =~ s/\s*<br\s*\/?>\s*/ /g;
$html =~ s/\s+/ /g;
($html, my @tableList) = extractTag($html, 'table');
foreach my $table (@tableList) {
	my $table = unwrapTag($table, 'table');
	$table = unwrapTag($table, 'thead');
	$table = unwrapTag($table, 'tbody');
	$table =~ s/^\s+//;
	$table =~ s/\s+$//;
	($table, my @trList) = extractTag($table, 'tr');
	print STDERR $table, "\n" if($table !~ /^\s*$/);
	foreach my $tr (@trList) {
		my $tr = unwrapTag($tr, 'tr');
		($tr, my @thList) = extractTag($tr, 'th');
		($tr, my @tdList) = extractTag($tr, 'td');
		$tr =~ s/^\s+//;
		$tr =~ s/\s+$//;
		print STDERR $tr, "\n" if($tr !~ /^\s*$/);
		if(@thList) {
			@thList = map {unwrapTag($_, 'th')} @thList;
			s/^\s+// foreach(@thList);
			s/\s+$// foreach(@thList);
			print join("\t", @thList), "\n";
		}
		if(@tdList) {
			@tdList = map {unwrapTag($_, 'td')} @tdList;
			s/^\s+// foreach(@tdList);
			s/\s+$// foreach(@tdList);
			print join("\t", @tdList), "\n";
		}
	}
}

sub unwrapTag {
	my ($html, $tag) = @_;
	while($html =~ /(<\/$tag>)/) {
		my ($closingTagStartIndex, $closingTagEndIndex) = ($-[1], $+[1]);
		if(substr($html, 0, $closingTagStartIndex) =~ /.*(<$tag)(>|\s)/s) {
			my ($openingTagStartIndex, $openingTagEndIndex) = ($-[1], $+[1]);
			$openingTagEndIndex += $+[1] if(substr($html, $openingTagEndIndex) =~ /^(\s+[^\s>]+\s*=\s*([^\s>]*|"[^"]*"|'[^']*')\s*)/);
			$openingTagEndIndex += $+[1] if(substr($html, $openingTagEndIndex) =~ /(>)/);
			my $closingTag = substr($html, $closingTagStartIndex, $closingTagEndIndex - $closingTagStartIndex, '');
			my $openingTag = substr($html, $openingTagStartIndex, $openingTagEndIndex - $openingTagStartIndex, '');
#			print STDERR "$openingTag\n";
#			print STDERR "$closingTag\n";
		} else {
			die 'error';
		}
	}
	return $html;
}

sub extractTag {
	my ($html, $tag) = @_;
	my @htmlList = ();
	while($html =~ /(<\/$tag>)/) {
		my ($closingTagStartIndex, $closingTagEndIndex) = ($-[1], $+[1]);
		if(substr($html, 0, $closingTagStartIndex) =~ /.*(<$tag)(>|\s)/s) {
			my ($openingTagStartIndex, $openingTagEndIndex) = ($-[1], $+[1]);
			$openingTagEndIndex += $+[1] if(substr($html, $openingTagEndIndex) =~ /^(\s+[^\s>]+\s*=\s*([^\s>]*|"[^"]*"|'[^']*')\s*)/);
			$openingTagEndIndex += $+[1] if(substr($html, $openingTagEndIndex) =~ /(>)/);
			push(@htmlList, substr($html, $openingTagStartIndex, $closingTagEndIndex - $openingTagStartIndex, ''));
		} else {
			die 'error';
		}
	}
	return ($html, @htmlList);
}
