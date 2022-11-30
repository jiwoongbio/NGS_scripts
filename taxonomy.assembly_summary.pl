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
	't' => \(my $replaceTaxid = ''),
	'n=i' => \(my $number = 0),
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

my %taxonIdHash = ();
my %taxonIdNumberHash = ();
my ($inputFile, $assemblySummaryFile) = @ARGV;
{
	open(my $reader, ($inputFile =~ /\.gz$/ ? "gzip -dc $inputFile |" : $inputFile)) or die "Can't open '$inputFile': $!";
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my $taxonId = $tokenList[$index];
		my @descendentTaxonIdList = ();
		if(my $taxon = $db->get_taxon(-taxonid => $taxonId)) {
			push(@descendentTaxonIdList, sort {$a <=> $b} map {$_->id} $db->get_all_Descendents($taxon));
		}
		$taxonIdHash{$_}->{$taxonId} = 1 foreach($taxonId, @descendentTaxonIdList);
		$taxonIdNumberHash{$taxonId} = 0;
	}
	close($reader);
}
{
	open(my $reader, ($assemblySummaryFile =~ /\.gz$/ ? "gzip -dc $assemblySummaryFile |" : $assemblySummaryFile)) or die "Can't open '$assemblySummaryFile': $!";
	my $line;
	chomp($line = <$reader>);
	print $line, "\n";
	chomp($line = <$reader>);
	print $line, "\n";
	$line =~ s/^# //;
	my @columnList = split(/\t/, $line, -1);
	while($line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		my @taxonIdList = ();
		foreach my $taxonId (sort {$a <=> $b} keys %{$taxonIdHash{$tokenHash{'taxid'}}}) {
			if($taxonIdNumberHash{$taxonId} < $number || $number == 0) {
				push(@taxonIdList, $taxonId);
				$taxonIdNumberHash{$taxonId} += 1;
			}
		}
		if(@taxonIdList) {
			$tokenHash{'taxid'} = join(',', @taxonIdList) if($replaceTaxid);
			print join("\t", @tokenHash{@columnList}), "\n";
		}
	}
	close($reader);
}
