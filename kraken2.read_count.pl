#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Taxonomy;
use List::Util qw(sum);
use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $directory = `readlink -f $0`);
$directory =~ s/\/[^\/]*$//;

my @excludeTaxonIdList = ();
GetOptions(
	'd=s' => \$directory,
	'r' => \(my $redownload = ''),
	't=s' => \(my $targetRank = ''),
	'T=s' => \(my $topRank = ''),
	'm=s' => \(my $moleculeNameDelimiter = ''),
	'e=i' => \@excludeTaxonIdList,
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
my %excludeTaxonIdHash = ();
foreach my $taxonId (@excludeTaxonIdList) {
	if(defined(my $taxon = $db->get_taxon(-taxonid => $taxonId))) {
		$excludeTaxonIdHash{$_} = 1 foreach($taxonId, map {$_->id} $db->get_all_Descendents($taxon));
	} else {
		print STDERR "Can't find taxon for '$taxonId'\n";
	}
}

my (@kraken2FileList) = @ARGV;
my %taxonIdCountHash = ();
{
	my %taxonIdTaxonIdHash = ();
	my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k3,3nr -k2,2n | uniq");
	foreach my $kraken2File (@kraken2FileList) {
		open(my $reader, ($kraken2File =~ /\.gz$/ ? "gzip -dc $kraken2File |" : $kraken2File)) or die "Can't open '$kraken2File': $!";
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line, -1);
			next if($tokenList[0] eq 'U');
			if($tokenList[2] =~ / \(taxid ([0-9]+)\)$/) {
				my $taxonId = $1;
				unless(defined($taxonIdTaxonIdHash{$taxonId})) {
					if(defined(my $taxon = $db->get_taxon(-taxonid => $taxonId))) {
						$taxonIdTaxonIdHash{$taxonId}->{$_->id} = 1 foreach(getTaxonList($taxon));
					} else {
						print STDERR "Can't find taxon for '$taxonId'\n";
					}
				}
				my %taxonIdKmerCountHash = ();
				$taxonIdKmerCountHash{$_->[0]} += $_->[1] foreach(map {[split(/:/, $_, 2)]} split(/ /, $tokenList[4]));
				next if(grep {$excludeTaxonIdHash{$_}} keys %taxonIdKmerCountHash);
				my $kmerCount = sum(0, @taxonIdKmerCountHash{grep {$taxonIdTaxonIdHash{$taxonId}->{$_}} keys %taxonIdKmerCountHash});
				my $moleculeName = getMoleculeName($tokenList[1]);
				print $writer join("\t", $moleculeName, $taxonId, $kmerCount), "\n" if($kmerCount > 0);
			} else {
				print STDERR "Can't parse '$tokenList[2]'\n";
			}
		}
		close($reader);
	}
	close($writer);
	my ($previousMoleculeName, $previousKmerCount) = ('');
	while(my $line = <$reader>) {
		chomp($line);
		my ($moleculeName, $taxonId, $kmerCount) = split(/\t/, $line, -1);
		if($previousMoleculeName ne $moleculeName) {
			($previousMoleculeName, $previousKmerCount) = ($moleculeName, $kmerCount);
		}
		$taxonIdCountHash{$taxonId} += 1 if($previousKmerCount == $kmerCount);
	}
	close($reader);
	waitpid($pid, 0);
}
my %taxonIdTaxonHash = ();
foreach my $taxonId (sort {$a <=> $b} keys %taxonIdCountHash) {
	my $taxon = $db->get_taxon(-taxonid => $taxonId);
	if($targetRank ne '') {
		$taxonIdTaxonHash{$_->id} = $_ foreach(grep {$_->rank eq $targetRank} getTaxonList($taxon));
	} else {
		$taxonIdTaxonHash{$taxon->id} = $taxon;
	}
}
foreach my $taxonId (sort {$a <=> $b} keys %taxonIdTaxonHash) {
	my $taxon = $taxonIdTaxonHash{$taxonId};
	my $count = 0;
	$count += $_ foreach(grep {defined} @taxonIdCountHash{map {$_->id} getTaxonList($taxon)});
	if($targetRank ne '') {
		$count += $_ foreach(grep {defined} @taxonIdCountHash{map {$_->id} $db->get_all_Descendents($taxon)});
	}
	print join("\t", $taxonId, $taxon->scientific_name, $taxon->rank, $count), "\n" if($count > 0);
}

sub getMoleculeName {
	my ($readName) = @_;
	my $moleculeName = $readName;
	$moleculeName =~ s/$moleculeNameDelimiter.*$// if($moleculeNameDelimiter ne '');
	return $moleculeName;
}

sub getTaxonList {
	my ($taxon) = @_;
	my @taxonList = ();
	while(defined($taxon)) {
		push(@taxonList, $taxon);
		return @taxonList if($topRank ne '' && $taxon->rank eq $topRank);
		$taxon = $taxon->ancestor;
	}
	@taxonList = () if($topRank ne '');
	return @taxonList;
}
