#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(max);
use Getopt::Long qw(:config no_ignore_case);

my @fastqFileList1 = ();
my @fastqFileList2 = ();
my @pairedEndLibraryList = ();
my @fastqFileListU = ();
my @singleEndLibraryList = ();
GetOptions(
	'h' => \(my $help = ''),
	'1=s' => \@fastqFileList1,
	'2=s' => \@fastqFileList2,
	'l=s' => \@pairedEndLibraryList,
	'U=s' => \@fastqFileListU,
	'L=s' => \@singleEndLibraryList,
	's=s' => \(my $sample = ''),
	'p=i' => \(my $threads = 8),
	'm=i' => \(my $memory = 32),
);
if($help || (scalar(@fastqFileList1) == 0 && scalar(@fastqFileList2) == 0 && scalar(@fastqFileListU) == 0)) {
	die <<EOF;

Usage:   perl chromosome_plasmid.assembly.mapping.pl {-1 paired.1.fastq -2 paired.2.fastq | -S single.fastq} [options]

Options: -h       display this help message
         -1 FILE  paired FASTQ file 1
         -2 FILE  paired FASTQ file 2
         -l STR   library name of paired FASTQ files
         -U FILE  single FASTQ file
         -L STR   library name of single FASTQ file
         -s STR   sample name [name of current directory]
         -p INT   number of threads [$threads]
         -m INT   memory limit in Gb [$memory]

EOF
}
my %libraryNumberHash = ();
my @libraryNumberList = ();
my @fastqFileListList = ();
if(@pairedEndLibraryList) {
	exit if(scalar(@pairedEndLibraryList) != scalar(@fastqFileList1) || scalar(@pairedEndLibraryList) != scalar(@fastqFileList2));
	foreach my $index (0 .. $#pairedEndLibraryList) {
		my $library = $pairedEndLibraryList[$index];
		$libraryNumberHash{$library} = scalar(keys %libraryNumberHash) + 1 unless($libraryNumberHash{$library});
		push(@libraryNumberList, $libraryNumberHash{$library});
		push(@fastqFileListList, [$fastqFileList1[$index], $fastqFileList2[$index]]);
	}
} else {
	exit if(scalar(@fastqFileList1) != scalar(@fastqFileList2));
	foreach my $index (0 .. max($#fastqFileList1, $#fastqFileList2)) {
		push(@fastqFileListList, [$fastqFileList1[$index], $fastqFileList2[$index]]);
	}
}
if(@singleEndLibraryList) {
	exit if(scalar(@singleEndLibraryList) != scalar(@fastqFileListU));
	foreach my $index (0 .. $#singleEndLibraryList) {
		my $library = $singleEndLibraryList[$index];
		$libraryNumberHash{$library} = scalar(keys %libraryNumberHash) + 1 unless($libraryNumberHash{$library});
		push(@libraryNumberList, $libraryNumberHash{$library});
		push(@fastqFileListList, [$fastqFileListU[$index]]);
	}
} else {
	foreach my $index (0 .. $#fastqFileListU) {
		push(@fastqFileListList, [$fastqFileListU[$index]]);
	}
}
if(@libraryNumberList) {
	exit if(scalar(@libraryNumberList) != scalar(@fastqFileListList));
}
if($sample eq '') {
	chomp(my $sample = `pwd`);
	$sample =~ s/^.*\///;
}

my @spadesInputList = ();
if(@libraryNumberList) {
	foreach my $index (0 .. $#fastqFileListList) {
		if(scalar(@{$fastqFileListList[$index]}) == 2) {
			push(@spadesInputList, "--pe$libraryNumberList[$index]-1 $fastqFileListList[$index]->[0]");
			push(@spadesInputList, "--pe$libraryNumberList[$index]-2 $fastqFileListList[$index]->[1]");
		} elsif(scalar(@{$fastqFileListList[$index]}) == 1) {
			push(@spadesInputList, "--s$libraryNumberList[$index] $fastqFileListList[$index]->[0]");
		}
	}
} else {
	foreach my $index (0 .. $#fastqFileListList) {
		if(scalar(@{$fastqFileListList[$index]}) == 2) {
			push(@spadesInputList, "-1 $fastqFileListList[$index]->[0]");
			push(@spadesInputList, "-2 $fastqFileListList[$index]->[1]");
		} elsif(scalar(@{$fastqFileListList[$index]}) == 1) {
			push(@spadesInputList, "-s $fastqFileListList[$index]->[0]");
		}
	}
}
system("rm -rf $sample.plasmidSPAdes");
system("time (spades.py --plasmid @spadesInputList -o $sample.plasmidSPAdes --threads $threads --memory $memory)");
my $noPlasmid = `grep 'No putative plasmid contigs found!' $sample.plasmidSPAdes/spades.log`;
if($noPlasmid eq '') {
	$noPlasmid = 1 if(-z "$sample.plasmidSPAdes/contigs.fasta");
}
if($noPlasmid eq '') {
	unless(-s "$sample.plasmidSPAdes/scaffolds.fasta") {
		print STDERR "$sample.plasmidSPAdes/scaffolds.fasta is not available.\n";
		$noPlasmid = 1;
	}
}
if($noPlasmid eq '') {
	system("time (sed 's/^>/>plasmid_/' $sample.plasmidSPAdes/scaffolds.fasta > $sample.plasmidSPAdes/plasmid_scaffolds.fasta)");
	system("time (bwa index $sample.plasmidSPAdes/plasmid_scaffolds.fasta)");
	foreach my $index (0 .. $#fastqFileListList) {
		system("time (bwa mem -t $threads -Y $sample.plasmidSPAdes/plasmid_scaffolds.fasta @{$fastqFileListList[$index]} | gzip > $sample.plasmidSPAdes/$index.sam.gz)");
		if(scalar(@{$fastqFileListList[$index]}) == 2) {
			system("time (samtools fastq -1 $sample.plasmidSPAdes/$index.unmapped.1.fastq -2 $sample.plasmidSPAdes/$index.unmapped.2.fastq -F 2306 -n $sample.plasmidSPAdes/$index.sam.gz)");
		} elsif(scalar(@{$fastqFileListList[$index]}) == 1) {
			system("time (samtools fastq -s $sample.plasmidSPAdes/$index.unmapped.fastq -F 2306 -n $sample.plasmidSPAdes/$index.sam.gz)");
		}
	}
	@spadesInputList = ();
	if(@libraryNumberList) {
		foreach my $index (0 .. $#fastqFileListList) {
			if(scalar(@{$fastqFileListList[$index]}) == 2) {
				push(@spadesInputList, "--pe$libraryNumberList[$index]-1 $sample.plasmidSPAdes/$index.unmapped.1.fastq");
				push(@spadesInputList, "--pe$libraryNumberList[$index]-2 $sample.plasmidSPAdes/$index.unmapped.2.fastq");
			} elsif(scalar(@{$fastqFileListList[$index]}) == 1) {
				push(@spadesInputList, "--s$libraryNumberList[$index] $sample.plasmidSPAdes/$index.unmapped.fastq");
			}
		}
	} else {
		foreach my $index (0 .. $#fastqFileListList) {
			if(scalar(@{$fastqFileListList[$index]}) == 2) {
				push(@spadesInputList, "-1 $sample.plasmidSPAdes/$index.unmapped.1.fastq");
				push(@spadesInputList, "-2 $sample.plasmidSPAdes/$index.unmapped.2.fastq");
			} elsif(scalar(@{$fastqFileListList[$index]}) == 1) {
				push(@spadesInputList, "-s $sample.plasmidSPAdes/$index.unmapped.fastq");
			}
		}
	}
}
system("rm -rf $sample.chromosomeSPAdes");
system("time (spades.py @spadesInputList -o $sample.chromosomeSPAdes --threads $threads --memory $memory)");
system("time (sed 's/^>/>chromosome_/' $sample.chromosomeSPAdes/scaffolds.fasta > $sample.chromosomeSPAdes/chromosome_scaffolds.fasta)");
system("time (bwa index $sample.chromosomeSPAdes/chromosome_scaffolds.fasta)");
foreach my $index (0 .. $#fastqFileListList) {
	system("time (bwa mem -t $threads -Y $sample.chromosomeSPAdes/chromosome_scaffolds.fasta @{$fastqFileListList[$index]} | gzip > $sample.chromosomeSPAdes/$index.sam.gz)");
}
if($noPlasmid eq '') {
	system("cat $sample.chromosomeSPAdes/chromosome_scaffolds.fasta $sample.plasmidSPAdes/plasmid_scaffolds.fasta > $sample.chromosome_plasmid.fasta");
	my @samFileList = ();
	foreach my $index (0 .. $#fastqFileListList) {
		system("time (remocon.pl -a -s $sample.chromosomeSPAdes/$index.sam.gz $sample.plasmidSPAdes/$index.sam.gz 2> $sample.chromosomeSPAdes/$index.remocon.log | gzip > $sample.chromosomeSPAdes/$index.remocon.sam.gz)");
		push(@samFileList, "$sample.chromosomeSPAdes/$index.remocon.sam.gz");
	}
	foreach my $index (0 .. $#fastqFileListList) {
		system("time (remocon.pl -a $sample.plasmidSPAdes/$index.sam.gz $sample.chromosomeSPAdes/$index.sam.gz 2> $sample.plasmidSPAdes/$index.remocon.log | gzip > $sample.plasmidSPAdes/$index.remocon.sam.gz)");
		push(@samFileList, "$sample.plasmidSPAdes/$index.remocon.sam.gz");
	}
	open(my $writer, "| gzip > $sample.chromosome_plasmid.sam.gz");
	concatenateSamFiles($writer, @samFileList);
	close($writer);
} else {
	system("cat $sample.chromosomeSPAdes/chromosome_scaffolds.fasta > $sample.chromosome_plasmid.fasta");
	my @samFileList = ();
	foreach my $index (0 .. $#fastqFileListList) {
		push(@samFileList, "$sample.chromosomeSPAdes/$index.sam.gz");
	}
	open(my $writer, "| gzip > $sample.chromosome_plasmid.sam.gz");
	concatenateSamFiles($writer, @samFileList);
	close($writer);
}
if($noPlasmid eq '') {
	system("rm $sample.plasmidSPAdes/plasmid_scaffolds.fasta*");
	foreach my $index (0 .. $#fastqFileListList) {
		system("rm $sample.plasmidSPAdes/$index.sam.gz");
		if(scalar(@{$fastqFileListList[$index]}) == 2) {
			system("rm $sample.plasmidSPAdes/$index.unmapped.1.fastq");
			system("rm $sample.plasmidSPAdes/$index.unmapped.2.fastq");
		} elsif(scalar(@{$fastqFileListList[$index]}) == 1) {
			system("rm $sample.plasmidSPAdes/$index.unmapped.fastq");
		}
		system("rm $sample.plasmidSPAdes/$index.remocon.sam.gz");
		system("rm $sample.plasmidSPAdes/$index.remocon.log");
	}
	system("rm $sample.chromosomeSPAdes/chromosome_scaffolds.fasta*");
	foreach my $index (0 .. $#fastqFileListList) {
		system("rm $sample.chromosomeSPAdes/$index.sam.gz");
		system("rm $sample.chromosomeSPAdes/$index.remocon.sam.gz");
		system("rm $sample.chromosomeSPAdes/$index.remocon.log");
	}
} else {
	system("rm $sample.chromosomeSPAdes/chromosome_scaffolds.fasta*");
	foreach my $index (0 .. $#fastqFileListList) {
		system("rm $sample.chromosomeSPAdes/$index.sam.gz");
	}
}

sub concatenateSamFiles {
	my ($writer, @samFileList) = @_;
	my @indexList = (0 .. $#samFileList);
	my @readerList = ();
	my @lineList = ();
	foreach my $index (@indexList) {
		open(my $reader, ($samFileList[$index] =~ /\.gz$/ ? "gzip -dc $samFileList[$index] |" : $samFileList[$index]));
		$readerList[$index] = $reader;
	}
	my @tagList = ();
	my %tagLineNumberHash = ();
	foreach my $index (@indexList) {
		my $reader = $readerList[$index];
		while($lineList[$index] = <$reader>) {
			chomp($lineList[$index]);
			if($lineList[$index] =~ /^(@\S*)/) {
				my $tag = $1;
				push(@tagList, $tag) unless(defined($tagLineNumberHash{$tag}));
				$tagLineNumberHash{$tag}->{$lineList[$index]} = scalar(keys %{$tagLineNumberHash{$tag}});
			} else {
				last;
			}
		}
	}
	foreach my $tag (@tagList) {
		my %lineNumberHash = %{$tagLineNumberHash{$tag}};
		print $writer "$_\n" foreach(sort {$lineNumberHash{$a} <=> $lineNumberHash{$b}} keys %lineNumberHash);
	}
	foreach my $index (@indexList) {
		my $reader = $readerList[$index];
		print $writer "$lineList[$index]\n";
		while($lineList[$index] = <$reader>) {
			chomp($lineList[$index]);
			print $writer "$lineList[$index]\n";
		}
	}
	foreach my $index (@indexList) {
		my $reader = $readerList[$index];
		close($reader);
	}
}
