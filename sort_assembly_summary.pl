#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Bio::DB::Taxonomy;
use Getopt::Long qw(:config no_ignore_case);

my @refseqCategoryList = ('reference genome', 'representative genome', 'na');
my @assemblyLevelList = ('Complete Genome', 'Chromosome', 'Scaffold', 'Contig');

my @assemblySummaryFilePathList = (
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt',
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt',
);

my @targetTaxonomyIdList = ();

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	'T=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'D=s' => \(my $taxonomyDirectory = '.'),
	'd' => \(my $redownload = ''),
	't=s' => \@targetTaxonomyIdList,
	'r=s' => \(my $ranks = 'species'),
	'g=f' => \(my $genomeAlignmentLengthRatioCutoff = 0),
	'n=i' => \(my $maximumNumber = 0),
);
if($help) {
	die <<EOF;

Usage:   perl sort_assembly_summary.pl [options] [assembly_summary.txt] [...]

Options: -h       display this help message
         -T DIR   temporary directory [$temporaryDirectory]
         -p INT   number of threads [$threads]
         -D DIR   taxonomy directory [$taxonomyDirectory]
         -d       redownload taxonomy data
         -t STR   comma-separated target NCBI taxonomy IDs or file
         -r STR   target taxonomy ranks [$ranks]
         -g FLOAT genome alignment length ratio cutoff of redundant assembly
         -n INT   maximum number per taxonomy [$maximumNumber]

EOF
}
my $processTemporaryDirectory = "$temporaryDirectory/$hostname.$$";
system("rm -rf $processTemporaryDirectory");
system("mkdir -p $processTemporaryDirectory");

if(not -r "$taxonomyDirectory/nodes.dmp" or not -r "$taxonomyDirectory/names.dmp" or $redownload) {
	my $URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
	my $file = "$taxonomyDirectory/taxdump.tar.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	system("cd $taxonomyDirectory; tar -zxf taxdump.tar.gz nodes.dmp");
	system("cd $taxonomyDirectory; tar -zxf taxdump.tar.gz names.dmp");
	system("rm -f $taxonomyDirectory/$_") foreach('nodes', 'parents', 'names2id', 'id2names');
}
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $taxonomyDirectory, -nodesfile => "$taxonomyDirectory/nodes.dmp", -namesfile => "$taxonomyDirectory/names.dmp");

my %rankHash = map {$_ => 1} split(/,/, $ranks);

my %targetTaxonomyIdHash = ();
foreach my $targetTaxonomyId (map {split(/,/, $_)} @targetTaxonomyIdList) {
	if(-r $targetTaxonomyId) {
		chomp(my @targetTaxonomyIdList = `cat $targetTaxonomyId`);
		$targetTaxonomyIdHash{$_} = 1 foreach(@targetTaxonomyIdList);
	} else {
		$targetTaxonomyIdHash{$targetTaxonomyId} = 1;
	}
}
@targetTaxonomyIdList = sort {$a <=> $b} keys %targetTaxonomyIdHash;
my %taxonomyIdTargetTaxonomyIdHash = ();
foreach my $targetTaxonomyId (@targetTaxonomyIdList) {
	if(my $taxon = $db->get_taxon(-taxonid => $targetTaxonomyId)) {
		$taxonomyIdTargetTaxonomyIdHash{$_}->{$targetTaxonomyId} = 1 foreach($targetTaxonomyId, map {$_->id} $db->get_all_Descendents($taxon));
	}
}

my %refseqCategoryIndexHash = map {$refseqCategoryList[$_] => $_} 0 .. $#refseqCategoryList;
my %assemblyLevelIndexHash = map {$assemblyLevelList[$_] => $_} 0 .. $#assemblyLevelList;

my (@assemblySummaryFileList) = @ARGV;
unless(@assemblySummaryFileList) {
	foreach my $assemblySummaryFilePath (@assemblySummaryFilePathList) {
		(my $assemblySummaryFile = $assemblySummaryFilePath) =~ s/^.*\///;
		$assemblySummaryFile = "$processTemporaryDirectory/$assemblySummaryFile";
		system("wget --no-verbose -O $assemblySummaryFile $assemblySummaryFilePath");
		push(@assemblySummaryFileList, $assemblySummaryFile);
	}
}
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6nr | cut -f1,7-");
my @headerLineList = ();
my @outputColumnList = ();
foreach my $assemblySummaryFileIndex (0 .. $#assemblySummaryFileList) {
	my $assemblySummaryFile = $assemblySummaryFileList[$assemblySummaryFileIndex];
	open(my $reader, ($assemblySummaryFile =~ /\.gz$/ ? "gzip -dc $assemblySummaryFile |" : $assemblySummaryFile)) or die "Can't open '$assemblySummaryFile': $!";
	my $line;
	chomp($line = <$reader>);
	push(@headerLineList, $line) if($assemblySummaryFileIndex == 0);
	chomp($line = <$reader>);
	push(@headerLineList, $line) if($assemblySummaryFileIndex == 0);
	$line =~ s/^# ?//;
	my @columnList = split(/\t/, $line, -1);
	@outputColumnList = @columnList if($assemblySummaryFileIndex == 0);
	while($line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		next if($tokenHash{'ftp_path'} eq 'na');
		my $taxonomyId = getTaxonomyId($tokenHash{'taxid'});
		next if($taxonomyId eq '');
		next if(%taxonomyIdTargetTaxonomyIdHash && !defined($taxonomyIdTargetTaxonomyIdHash{$taxonomyId}));
		my $refseqCategoryIndex = $refseqCategoryIndexHash{$tokenHash{'refseq_category'}};
		my $assemblyLevelIndex = $assemblyLevelIndexHash{$tokenHash{'assembly_level'}};
		(my $assemblyAccessionNumber = $tokenHash{'assembly_accession'}) =~ s/^[^0-9]*//;
		my $assemblyAccessionNumberVersion = 0;
		if($assemblyAccessionNumber =~ s/\.([0-9]+)$//) {
			$assemblyAccessionNumberVersion = $1;
		}
		print $writer join("\t", $taxonomyId, $refseqCategoryIndex, $assemblyLevelIndex, $assemblySummaryFileIndex, $assemblyAccessionNumber, $assemblyAccessionNumberVersion, @tokenHash{@outputColumnList}), "\n";
	}
	close($reader);
}
close($writer);
{
	foreach my $line (@headerLineList) {
		print $line, "\n";
	}
	my $previousTaxonomyId = 0;
	my $count = 0;
	my %assemblyAccessionHash = ();
	my @assemblyPrefixList = ();
	my %childPidHash = ();
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(my $taxonomyId, @tokenHash{@outputColumnList}) = split(/\t/, $line, -1);
		if($previousTaxonomyId != $taxonomyId) {
			$previousTaxonomyId = $taxonomyId;
			$count = 0;
			%assemblyAccessionHash = ();
			@assemblyPrefixList = ();
		}
		if($maximumNumber == 0 || $count < $maximumNumber) {
			unless($assemblyAccessionHash{$tokenHash{'assembly_accession'}}) {
				unless(isRedundantAssembly(%tokenHash)) {
					print join("\t", @tokenHash{@outputColumnList}), "\n";
					$count += 1;
				}
			}
		}
		$assemblyAccessionHash{$tokenHash{'assembly_accession'}} = 1;
		$assemblyAccessionHash{$tokenHash{'gbrs_paired_asm'}} = 1 if($tokenHash{'paired_asm_comp'} eq 'identical');
	}
	sub isRedundantAssembly {
		my %tokenHash = @_;
		if($genomeAlignmentLengthRatioCutoff > 0) {
			(my $assemblyPrefix = $tokenHash{'ftp_path'}) =~ s/^.*\///;
			my $genomeLength = 0;
			for(my $try = 0; $try < 5; $try += 1) {
				system("wget --no-verbose -O - $tokenHash{'ftp_path'}/$assemblyPrefix\_genomic.fna.gz | gzip -d > $processTemporaryDirectory/$assemblyPrefix.fasta");
				open(my $reader, "$processTemporaryDirectory/$assemblyPrefix.fasta");
				while(my $line = <$reader>) {
					chomp($line);
					$genomeLength += length($line) unless($line =~ /^>/);
				}
				close($reader);
				last if($genomeLength > 0);
			}
			return 1 if($genomeLength == 0);
			foreach my $previousAssemblyPrefix (@assemblyPrefixList) {
				if(my $childPid = fork()) {
					$childPidHash{$childPid} = 1;
				} else {
					my $prefix = "$processTemporaryDirectory/$$";
					system("nucmer --mum --prefix=$prefix $processTemporaryDirectory/$assemblyPrefix.fasta $processTemporaryDirectory/$previousAssemblyPrefix.fasta");
					system("delta-filter -q -r $prefix.delta > $prefix.filter.delta");
					my $genomeAlignmentLength = 0;
					open(my $reader, "show-coords -H -r -T $prefix.filter.delta |");
					while(my $line = <$reader>) {
						chomp($line);
						my ($start1, $end1, $start2, $end2, $length1, $length2) = split(/\t/, $line, -1);
						$genomeAlignmentLength += $length1;
					}
					close($reader);
					system("rm $prefix.delta $prefix.filter.delta");
					my $genomeAlignmentLengthRatio = $genomeAlignmentLength / $genomeLength;
					open(my $writer, "> $prefix.genomeAlignmentLengthRatio.txt");
					print $writer join("\t", $assemblyPrefix, $previousAssemblyPrefix, $genomeAlignmentLengthRatio), "\n";
					close();
					exit(0);
				}
				if(isRedundantAssembly_wait($threads)) {
					return 1;
				}
			}
			if(isRedundantAssembly_wait()) {
				return 1;
			}
			push(@assemblyPrefixList, $assemblyPrefix);
		}
		return '';
	}
	sub isRedundantAssembly_wait {
		my ($number) = (@_, 1);
		my $isRedundantAssembly = '';
		while(scalar(keys %childPidHash) >= $number) {
			my $childPid = wait();
			if($childPidHash{$childPid}) {
				open(my $reader, "$processTemporaryDirectory/$childPid.genomeAlignmentLengthRatio.txt");
				while(my $line = <$reader>) {
					chomp($line);
					print STDERR $line, "\n";
					my ($assemblyPrefix, $previousAssemblyPrefix, $genomeAlignmentLengthRatio) = split(/\t/, $line, -1);
					if($genomeAlignmentLengthRatio > $genomeAlignmentLengthRatioCutoff) {
						$isRedundantAssembly = 1;
						$number = 1;
					}
				}
				close($reader);
				system("rm $processTemporaryDirectory/$childPid.genomeAlignmentLengthRatio.txt");
				delete $childPidHash{$childPid};
			}
		}
		return $isRedundantAssembly;
	}
}
close($reader);
waitpid($pid, 0);

system("rm -rf $processTemporaryDirectory");

sub getTaxonomyId {
	my ($taxonomyId) = @_;
	my $taxon = $db->get_taxon(-taxonid => $taxonomyId);
	while(defined($taxon)) {
		return $taxon->id if($rankHash{$taxon->rank});
		$taxon = $taxon->ancestor;
	}
	return '';
}
