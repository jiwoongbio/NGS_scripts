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
	'd=s' => \$directory,
	'r' => \(my $redownload = ''),
	't=s' => \(my $targetRank = ''),
	'T=s' => \(my $toRank = ''),
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

my (@kraken2FileList) = @ARGV;
my %taxonIdCountHash = ();
foreach my $kraken2File (@kraken2FileList) {
	open(my $reader, ($kraken2File =~ /\.gz$/ ? "gzip -dc $kraken2File |" : $kraken2File)) or die "Can't open '$kraken2File': $!";
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		if($tokenList[2] =~ / \(taxid ([0-9]+)\)$/) {
			$taxonIdCountHash{$1} += 1;
		} else {
			print STDERR "Can't parse '$tokenList[2]'\n";
		}
	}
	close($reader);
}
foreach my $taxonId (sort {$a <=> $b} keys %taxonIdCountHash) {
	if(defined(my $taxon = $db->get_taxon(-taxonid => $taxonId))) {
		next if($targetRank ne '' && $targetRank ne $taxon->rank);
		my $count = 0;
		$count += $_ foreach(grep {defined} @taxonIdCountHash{getTaxonIdList($taxonId)});
		if($targetRank ne '') {
			$count += $_ foreach(grep {defined} @taxonIdCountHash{map {$_->id} $db->get_all_Descendents($taxon)});
		}
		print join("\t", $taxonId, $taxon->scientific_name, $taxon->rank, $count), "\n";
	} else {
		print STDERR "Can't find taxon for '$taxonId'\n";
	}
}

sub getTaxonIdList {
	my ($taxonId) = @_;
	my $taxon = $db->get_taxon(-taxonid => $taxonId);
	my @taxonIdList = ();
	while(defined($taxon)) {
		push(@taxonIdList, $taxon->id);
		return @taxonIdList if($toRank ne '' && $toRank eq $taxon->rank);
		$taxon = $taxon->ancestor;
	}
	@taxonIdList = () if($toRank ne '');
	return @taxonIdList;
}
