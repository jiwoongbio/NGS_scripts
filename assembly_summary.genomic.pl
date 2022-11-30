#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my (@assemblySummaryFileList) = @ARGV;
my %accessionHash = ();
foreach my $assemblySummaryFile (@assemblySummaryFileList) {
	open(my $reader, ($assemblySummaryFile =~ /\.gz$/ ? "gzip -dc $assemblySummaryFile |" : $assemblySummaryFile)) or die "Can't open '$assemblySummaryFile': $!";
	my $line;
	chomp($line = <$reader>);
	chomp($line = <$reader>);
	$line =~ s/^# //;
	my @columnList = split(/\t/, $line, -1);
	while($line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		my $accession = $tokenHash{'assembly_accession'};
		next if($accessionHash{$accession});
		(my $prefix = $tokenHash{'ftp_path'}) =~ s/^.*\///;
		open(my $reader, "wget --no-verbose -O - $tokenHash{'ftp_path'}/${prefix}_genomic.fna.gz | gzip -d |");
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ s/^>//) {
				print ">$accession|$line\n";
			} else {
				print "$line\n";
			}
		}
		close($reader);
		$accessionHash{$accession} = 1;
	}
	close($reader);
}
