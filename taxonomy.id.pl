#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use URI::Escape;
use Bio::DB::Taxonomy;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $directory = `readlink -f $0`);
$directory =~ s/\/[^\/]*$//;
GetOptions(
	'i=i' => \(my $index = 0),
	'd=s' => \$directory,
	'r' => \(my $redownload = ''),
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

my ($inputFile) = @ARGV;
open(my $reader, ($inputFile =~ /\.gz$/ ? "gzip -dc $inputFile |" : $inputFile)) or die "Can't open '$inputFile': $!";
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my $taxonName = $tokenList[$index];
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
	push(@taxonIdList, '') unless(@taxonIdList);
	print join("\t", $line, $_), "\n" foreach(@taxonIdList);
}
close($reader);

sub esearch {
	my ($db, $term) = @_;
	my $encodedTerm = uri_escape($term);
	chomp(my @idList = `esearch.pl -e $db '$encodedTerm'`);
	return @idList;
}
