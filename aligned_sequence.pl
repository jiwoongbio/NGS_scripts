#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);
use Bio::DB::Fasta;

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	't=s' => \$temporaryDirectory,
);
my $temporaryPrefix = "$temporaryDirectory/.$hostname.$$";
my ($referenceFastaFile, $queryFastaFile) = @ARGV;
my $referenceDB = Bio::DB::Fasta->new($referenceFastaFile);
my $queryDB = Bio::DB::Fasta->new($queryFastaFile);
system("nucmer --mum --prefix=$temporaryPrefix $referenceFastaFile $queryFastaFile");
system("delta-filter -q -r $temporaryPrefix.delta > $temporaryPrefix.filter.delta");
system("show-coords -H -r -T $temporaryPrefix.filter.delta > $temporaryPrefix.coords");
{
	open(my $reader, "$temporaryPrefix.coords");
	my ($name, $sequence, $nextPosition) = ('', '', 1);
	while(my $line = <$reader>) {
		chomp($line);
		my ($referenceStart, $referenceEnd, $queryStart, $queryEnd, $referenceLength, $queryLength, $identity, $referenceName, $queryName) = split(/\t/, $line, -1);
		if($name ne $referenceName) {
			if($name ne '') {
				$sequence .= $referenceDB->seq($name, $nextPosition, $referenceDB->length($name)) if($nextPosition <= $referenceDB->length($name));
				print ">$name\n";
				print "$sequence\n";
			}
			($name, $sequence, $nextPosition) = ($referenceName, '', 1);
		}
		if($nextPosition <= $referenceStart) {
			$sequence .= $referenceDB->seq($name, $nextPosition, $referenceStart - 1) if($nextPosition <= $referenceStart - 1);
			$sequence .= $queryDB->seq($queryName, $queryStart, $queryEnd);
			$nextPosition = $referenceEnd + 1;
		}
	}
	if($name ne '') {
		$sequence .= $referenceDB->seq($name, $nextPosition, $referenceDB->length($name)) if($nextPosition <= $referenceDB->length($name));
		print ">$name\n";
		print "$sequence\n";
	}
	close($reader);
}
system("rm $temporaryPrefix.delta $temporaryPrefix.filter.delta $temporaryPrefix.coords");
