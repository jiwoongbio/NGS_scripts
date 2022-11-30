#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Taxonomy;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $directory = `readlink -f $0`);
$directory =~ s/\/[^\/]*$//;
GetOptions(
	'i=i' => \(my $index = 0),
	'd=s' => \$directory,
	'r' => \(my $redownload = ''),
	't=s' => \(my $taxonNameIdFile = ''),
	's' => \(my $isSecondColumnTaxonId = ''),
	'u' => \(my $allowUnmatchedTaxonomy = ''),
	'd' => \(my $allowDisorder = ''),
	'r=s' => \(my $rank = ''),
);
if(not -r "$directory/nodes.dmp" or not -r "$directory/names.dmp" or $redownload) {
	my $URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
	my $file = "$directory/taxdump.tar.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	system("cd $directory; tar -zxf taxdump.tar.gz nodes.dmp");
	system("cd $directory; tar -zxf taxdump.tar.gz names.dmp");
	system("rm -f $directory/$_") foreach('nodes', 'parents', 'names2id', 'id2names');
}
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $directory, -nodesfile => "$directory/nodes.dmp", -namesfile => "$directory/names.dmp");

my %taxonNameIdHash = ();
if($taxonNameIdFile ne '') {
	open(my $reader, $taxonNameIdFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my ($taxonName, $taxonIds) = @tokenList;
		my @taxonIdList = split(/,/, $taxonIds);
		$taxonNameIdHash{$taxonName}->{$_} = 1 foreach(@taxonIdList);
	}
	close($reader);
}
my %taxonNameIdListHash = ();
foreach my $taxonName (keys %taxonNameIdHash) {
	my @taxonIdList = sort {$a <=> $b} keys %{$taxonNameIdHash{$taxonName}};
	$taxonNameIdListHash{$taxonName} = \@taxonIdList;
}

my ($silvaFile) = @ARGV;
open(my $reader, $silvaFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my @taxonNameList = split(/[;|]/, $tokenList[0]);
	s/^.*__// foreach(@taxonNameList);
	s/_/ /g foreach(@taxonNameList);
	my @taxonIdListList = ();
	foreach my $taxonName (@taxonNameList) {
		if($taxonNameIdFile ne '') {
			$taxonNameIdListHash{$taxonName} = [] unless(defined($taxonNameIdListHash{$taxonName}));
			push(@taxonIdListList, $taxonNameIdListHash{$taxonName});
		} else {
			my @taxonIdList = $db->get_taxonids($taxonName);
			unless(@taxonIdList) {
				@taxonIdList = esearch('taxonomy', "\"$taxonName\"");
				unless(@taxonIdList) {
					@taxonIdList = esearch('taxonomy', join(' AND ', map {"\"$_\""} split(/\s+/, $taxonName))) if($taxonName =~ /\s+/);
				}
			}
			if(scalar(@taxonIdList) > 1 && scalar(my @filteredTaxonIdList = grep {lc($db->get_taxon(-taxonid => $_)->scientific_name) eq lc($taxonName)} @taxonIdList) == 1) {
				@taxonIdList = @filteredTaxonIdList;
			}
			push(@taxonIdListList, \@taxonIdList);
		}
	}
	my @scoreTaxonIdListList = ();
	foreach my $taxonId ($isSecondColumnTaxonId ? $tokenList[1] : @{$taxonIdListList[-1]}) {
		my @taxonIdList = ();
		if(defined(my $taxon = $db->get_taxon(-taxonid => $taxonId))) {
			if($rank eq '' || $taxon->rank eq $rank) {
				unshift(@taxonIdList, $taxon->id);
				unshift(@taxonIdList, $taxon->id) while(defined($taxon = $taxon->ancestor()));
			}
		}
		my %taxonIdIndexHash = map {$taxonIdList[$_] => $_} 0 .. $#taxonIdList;
		my @indexListList = ();
		foreach(@taxonIdListList) {
			my @indexList = grep {defined($_)} @taxonIdIndexHash{@$_};
			push(@indexList, '') if(scalar(@indexList) == 0);
			push(@indexListList, \@indexList);
		}
		@indexListList = getCombinationList(@indexListList);
		foreach(@indexListList) {
			my @indexList = @$_;
			next if($allowUnmatchedTaxonomy eq '' && grep {$_ eq ''} @indexList);
			if(my @matchedIndexList = grep {$_ ne ''} @indexList) {
				next if($allowDisorder eq '' && grep {$matchedIndexList[$_ - 1] > $matchedIndexList[$_]} 1 .. $#matchedIndexList);
				my $score = 0;
				$score += $_ foreach(map {$indexList[$_] ne '' ? 2 ** $_ : 0} 0 .. $#indexList);
				push(@scoreTaxonIdListList, [$score, map {$_ eq '' ? '' : $taxonIdList[$_]} @indexList]);
			}
		}
	}
	if(@scoreTaxonIdListList) {
		foreach(@scoreTaxonIdListList) {
			my ($score, @taxonIdList) = @$_;
			print join("\t", $line, $score, join('', map {"$_;"} @taxonIdList)), "\n";
		}
	} else {
		print join("\t", $line, '', ''), "\n";
	}
}
close($reader);

sub getCombinationList {
	my (@tokenListList) = @_;
	my @combinationList = ();
	if(my ($index) = grep {ref($tokenListList[$_])} 0 .. $#tokenListList) {
		foreach(@{$tokenListList[$index]}) {
			push(@combinationList, getCombinationList(@tokenListList[0 .. ($index - 1)], $_, @tokenListList[($index + 1) .. $#tokenListList]));
		}
	} else {
		push(@combinationList, \@tokenListList);
	}
	return @combinationList;
}

sub esearch {
	my ($db, $term) = @_;
	my $encodedTerm = uri_escape($term);
	chomp(my @idList = `esearch.pl -e $db '$encodedTerm'`);
	return @idList;
}
