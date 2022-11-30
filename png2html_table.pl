#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);
use GD;

my ($pngFile, $numberOfColumns) = @ARGV;
my $image = GD::Image->newFromPng($pngFile);

my ($height, $width) = $image->getBounds;
my $denominator = ($width + 1) / $numberOfColumns;

my $rgb;
foreach my $x (0 .. $width) {
	foreach my $y (0 .. $height) {
		my ($index) = $image->getPixel($x, $y);
		my ($r, $g, $b) = $image->rgb($index);
		my $col = int($x / $denominator);
		my $row = int($y / $denominator);
		push(@{$rgb->[$row]->[$col]->{'r'}}, $r);
		push(@{$rgb->[$row]->[$col]->{'g'}}, $g);
		push(@{$rgb->[$row]->[$col]->{'b'}}, $b);
	}
}

print "<table>\n";
foreach my $row (0 .. scalar(@$rgb) - 1) {
	my @tdList = ();
	foreach my $col (0 .. scalar(@{$rgb->[$row]}) - 1) {
		my $r = mean(@{$rgb->[$row]->[$col]->{'r'}});
		my $g = mean(@{$rgb->[$row]->[$col]->{'g'}});
		my $b = mean(@{$rgb->[$row]->[$col]->{'b'}});
		(my $color = sprintf('#%02x%02x%02x', int($r), int($g), int($b))) =~ tr/a-z/A-Z/;
		push(@tdList, "<td style=\"background-color: $color;\"></td>");
	}
	print join('', '<tr>', @tdList, '</tr>'), "\n";
}
print "</table>\n";

sub mean {
	my @valueList = @_;
	return sum(@valueList) / scalar(@valueList);
}
