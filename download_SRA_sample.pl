#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Time::HiRes;
use URI::Escape;
use XML::LibXML;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'e' => \(my $termEncoded = ''),
	'mindate=s' => \(my $mindate = ''),
	'maxdate=s' => \(my $maxdate = ''),
	'E' => \(my $downloadFromEBI = ''),
	'x' => \(my $experimentInsteadOfSample = ''),
	'm' => \(my $doNotMerge = ''),
	'F' => \(my $origfmt = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] SRA_search_term

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -e       term is encoded
         -mindate STR   limit search by minimum date
         -maxdate STR   limit search by maximum date
         -E       download from EBI
         -x       experiment instead of sample
         -m       do not merge
         -F       Defline contains only original sequence name

EOF
}
my $temporaryPrefix = "$temporaryDirectory/$hostname.$$";
system("rm -f $temporaryPrefix.*");
{
	my %pidHash = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			$pidHash{$pid} = 1;
		} else {
			open($writer, "> $temporaryPrefix.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(keys %pidHash) >= $number) {
			my $pid = wait();
			if($pidHash{$pid}) {
				open(my $reader, "$temporaryPrefix.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryPrefix.$pid");
				delete $pidHash{$pid};
			}
		}
	}
	sub forkPrint {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parentWriter)) {
			print $parentWriter @_;
		} else {
			print @_;
		}
	}
}
my ($term) = @ARGV;
my $db = 'sra';
my @idList = esearch($db, $term);
my %layoutSampleRunListHash = ();
if(@idList) {
	my $retmax = 100;
	my $rettype = 'full';
	my $retmode = 'xml';
	for(my $index = 0; $index < scalar(@idList); $index += $retmax) {
		my $id = join(',', grep {defined} @idList[$index .. $index + $retmax - 1]);
		my $xmlString = efetch($db, $id, $rettype, $retmode);
		my $dom = XML::LibXML->load_xml(string => $xmlString);
		my $root = $dom->documentElement();
		foreach my $experimentPackageNode (getChildNodeList($root, 'EXPERIMENT_PACKAGE')) {
			my @layoutList = map {$_->nodeName} map {$_->childNodes()} getChildNodeList($experimentPackageNode, 'EXPERIMENT', 'DESIGN', 'LIBRARY_DESCRIPTOR', 'LIBRARY_LAYOUT');
			my @sampleList = map {$_->getAttribute('accession')} getChildNodeList($experimentPackageNode, 'SAMPLE');
			@sampleList = map {$_->getAttribute('accession')} getChildNodeList($experimentPackageNode, 'EXPERIMENT') if($experimentInsteadOfSample);
			my @runList = map {$_->getAttribute('accession')} getChildNodeList($experimentPackageNode, 'RUN_SET', 'RUN');
			foreach my $layout (@layoutList) {
				foreach my $sample (@sampleList) {
					push(@{$layoutSampleRunListHash{$layout}->{$sample}}, @runList);
				}
			}
		}
	}
}
foreach my $layout (sort keys %layoutSampleRunListHash) {
	system("mkdir -p $layout");
	open(my $writer, "> $layout/sample.txt");
	foreach my $sample (sort keys %{$layoutSampleRunListHash{$layout}}) {
		print $writer "$sample\n";
		if($threads == 1) {
			download($layout, $sample);
		} else {
			forkPrintSubroutine(\&download, $layout, $sample);
		}
	}
	close($writer);
}
forkPrintWait();

sub download {
	my ($layout, $sample) = @_;
	system("mkdir -p $layout/$sample");
	my @runList = sort @{$layoutSampleRunListHash{$layout}->{$sample}};
	if($downloadFromEBI) {
		if($layout eq 'SINGLE') {
			foreach my $run (@runList) {
				system(sprintf("cd $layout/$sample; wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/$run/$run\_1.fastq.gz", substr($run, 0, 6))) unless(-s "$layout/$sample/$run\_1.fastq.gz");
			}
		}
		if($layout eq 'PAIRED') {
			foreach my $run (@runList) {
				system(sprintf("cd $layout/$sample; wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/$run/$run\_1.fastq.gz", substr($run, 0, 6))) unless(-s "$layout/$sample/$run\_1.fastq.gz");
				system(sprintf("cd $layout/$sample; wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/$run/$run\_2.fastq.gz", substr($run, 0, 6))) unless(-s "$layout/$sample/$run\_2.fastq.gz");
			}
		}
	} else {
		foreach my $run (@runList) {
			if($origfmt) {
				system("cd $layout/$sample; fastq-dump --gzip --split-files --log-level err --origfmt $run");
			} else {
				system("cd $layout/$sample; fastq-dump --gzip --split-files --log-level err $run");
			}
			system("rm -f ~/ncbi/public/sra/$run ~/ncbi/public/sra/$run.*");
		}
	}
	if($doNotMerge eq '') {
		if($layout eq 'SINGLE') {
			my @fileList = map {"$layout/$sample/$_\_1.fastq.gz"} @runList;
			system("for file in @fileList; do gzip -dc \$file; done | gzip > $layout/$sample/$sample.1.fastq.gz");
			system("rm @fileList");
		}
		if($layout eq 'PAIRED') {
			my @fileList = map {"$layout/$sample/$_\_1.fastq.gz"} @runList;
			system("for file in @fileList; do gzip -dc \$file; done | gzip > $layout/$sample/$sample.1.fastq.gz");
			system("rm @fileList");
		}
		if($layout eq 'PAIRED') {
			my @fileList = map {"$layout/$sample/$_\_2.fastq.gz"} @runList;
			system("for file in @fileList; do gzip -dc \$file; done | gzip > $layout/$sample/$sample.2.fastq.gz");
			system("rm @fileList");
		}
	}
}

sub efetch {
	my ($db, $id, $rettype, $retmode) = @_;
	my $output = `efetch.pl $db $id $rettype $retmode`;
	return $output;
}

sub esearch {
	my ($db, $term) = @_;
	my $encodedTerm = $termEncoded ? $term : uri_escape($term);
	chomp(my @idList = `esearch.pl -e -mindate '$mindate' -maxdate '$maxdate' $db '$encodedTerm'`);
	return @idList;
}

sub getChildNodeList {
	my ($node, @childNodeTagNameList) = @_;
	my @childNodeList = ();
	if(@childNodeTagNameList) {
		foreach my $childNode ($node->getChildrenByTagName(shift @childNodeTagNameList)) {
			push(@childNodeList, getChildNodeList($childNode, @childNodeTagNameList));
		}
	} else {
		push(@childNodeList, $node);
	}
	return @childNodeList;
}
