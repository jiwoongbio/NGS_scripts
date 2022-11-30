#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;

foreach my $fastaFile (@ARGV) {
	my $db = Bio::DB::Fasta->new($fastaFile);
}
