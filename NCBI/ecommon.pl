#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

chomp(my $directory = `readlink -f $0`);
$directory =~ s/\/[^\/]*$//;

use Time::HiRes;

unless(-e "$directory/tmp") {
	system("mkdir -p $directory/tmp");
	system("chmod 777 $directory/tmp");
}
my $lockFile = "$directory/tmp/eutils.lock";
sleep(1) while(-e $lockFile);
system("touch $lockFile") == 0 or die "Failed to lock\n";

if(-e (my $file = "$directory/tmp/eutils.api_key")) {
	open(my $reader, $file);
	chomp($ENV{'NCBI_API_KEY'} = <$reader>);
	close($reader);
}

my ($programQueryString) = @ARGV;
my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $apiKey = $ENV{'NCBI_API_KEY'};
my $credit = 100;
while($credit > 0) {
	$credit -= 1;
	my $exitStatus;
	if(defined($apiKey)) {
		$exitStatus = system("wget --no-verbose --no-check-certificate -O - '$baseURL/$programQueryString&api_key=$apiKey'");
		Time::HiRes::usleep(105);
	} else {
		$exitStatus = system("wget --no-verbose --no-check-certificate -O - '$baseURL/$programQueryString'");
		Time::HiRes::usleep(350);
	}
	last if($exitStatus == 0);
}
die if($credit == 0);

system("rm $lockFile") == 0 or die "Failed to unlock\n";
