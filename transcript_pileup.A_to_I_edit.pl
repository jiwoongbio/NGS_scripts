#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($transcriptPileupFile, $transcriptVariantFile) = @ARGV;
my %variantChromosomePositionHash = ();
if(defined($transcriptVariantFile)) {
	open(my $reader, ($transcriptVariantFile =~ /\.gz$/ ? "gzip -dc $transcriptVariantFile |" : $transcriptVariantFile));
	while(my $line = <$reader>) {
		chomp($line);
		my ($chromosome, $position, $refBase, $altBase) = split(/\t/, $line, -1);
		next if($refBase ne 'A');
		next if($altBase ne 'G');
		$variantChromosomePositionHash{$chromosome}->{$position} = 1;
	}
	close($reader);
}
my $countA = 0;
my $countG = 0;
open(my $reader, ($transcriptPileupFile =~ /\.gz$/ ? "gzip -dc $transcriptPileupFile |" : $transcriptPileupFile));
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split("\t", $line, -1);
	my ($chromosome, $position, $refBase) = @tokenList;
	$refBase =~ tr/a-z/A-Z/;
	next if($refBase ne 'A');
	my %readBaseCountHash = ();
	for(my $index = 3; $index < scalar(@tokenList); $index += 3) {
		my ($readDepth, $readBases, $baseQualities) = @tokenList[$index .. $index + 2];
		foreach my $readBase (getReadBaseList($readBases)) {
			$readBase =~ s/^\^.//;
			$readBase =~ s/\$$//;
			next if($readBase eq '>');
			next if($readBase eq '<');
			next if($readBase =~ /^\*/);
			$readBase =~ s/^[.,]/$refBase/;
			$readBase =~ tr/a-z/A-Z/;
			$readBaseCountHash{$readBase} += 1;
		}
	}
	$countA += $readBaseCountHash{'A'} if($readBaseCountHash{'A'});
	unless($variantChromosomePositionHash{$chromosome}->{$position}) {
		$countG += $readBaseCountHash{'G'} if($readBaseCountHash{'G'});
	}
}
close($reader);
print join("\t", $countA, $countG, $countG / $countA), "\n";

sub getReadBaseList {
	my ($readBases) = @_;
	my @readBaseList = ();
	while($readBases ne '') {
		my $readBase = '';
		$readBase .= $1 if($readBases =~ s/^(\^.)//);
		$readBases =~ s/^([.A-Z>,a-z<*])//;
		$readBase .= $1;
		while($readBases =~ s/^([+-])([0-9]+)//) {
			$readBase .= "$1$2";
			$readBase .= substr($readBases, 0, $2, '');
		}
		$readBase .= $1 if($readBases =~ s/^(\$)//);
		push(@readBaseList, $readBase);
	}
	return @readBaseList;
}
