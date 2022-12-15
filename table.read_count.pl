#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);

my $threads = 8;
my @columnList = ();
my %columnHash = ();
my %columnSampleCountHash = ();
my (@sampleList) = @ARGV;
chomp(@sampleList = `cat sample.txt`) if(scalar(@sampleList) == 0);
foreach my $sample (@sampleList) {
	chomp(my $directory = `readlink -f $sample`);
	(my $prefix = $directory) =~ s/^.*\///;
	my $directoryPrefix = "$directory/$prefix";
	{
		addColumn('Sample');
		$columnSampleCountHash{'Sample'}->{$sample} = $sample;
	}
	if(-s "$directoryPrefix.read_length.read_count.txt") {
		addColumn('Read length', 'Total reads');
		chomp($columnSampleCountHash{'Read length'}->{$sample} = `cut -f2 $directoryPrefix.read_length.read_count.txt | sort -nr | head -n1`);
		chomp($columnSampleCountHash{'Total reads'}->{$sample} = `cut -f3 $directoryPrefix.read_length.read_count.txt | awk '{a += \$o} END {print a}'`);
	}
	if(-s "$directoryPrefix.trimgalore.read_length.read_count.txt") {
		addColumn('QC-passed reads', '% QC-passed reads');
		chomp($columnSampleCountHash{'QC-passed reads'}->{$sample} = `cut -f3 $directoryPrefix.trimgalore.read_length.read_count.txt | awk '{a += \$o} END {print a}'`);
		$columnSampleCountHash{'% QC-passed reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'QC-passed reads'}->{$sample} / $columnSampleCountHash{'Total reads'}->{$sample} * 100);
	}
	if(-s "$directoryPrefix.trimmed.read_length.read_count.txt") {
		addColumn('QC-passed reads', '% QC-passed reads');
		chomp($columnSampleCountHash{'QC-passed reads'}->{$sample} = `cut -f3 $directoryPrefix.trimmed.read_length.read_count.txt | awk '{a += \$o} END {print a}'`);
		$columnSampleCountHash{'% QC-passed reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'QC-passed reads'}->{$sample} / $columnSampleCountHash{'Total reads'}->{$sample} * 100);
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.transcript.remocon.log", "$directoryPrefix.transcriptome.remocon.log")) {
		chomp(my @typeCountList = `cat $file`);
		my %typeCountHash = map {$_->[0] => $_->[1]} map {[split(/: /, $_, 2)]} @typeCountList;
		my $count = sum(values %typeCountHash);
		if($count * 2 == $columnSampleCountHash{'QC-passed reads'}->{$sample}) {
			$typeCountHash{$_} = $typeCountHash{$_} * 2 foreach(keys %typeCountHash);
		} elsif($count != $columnSampleCountHash{'QC-passed reads'}->{$sample}) {
			print STDERR join("\t", $sample, 'different read count', 'transcript'), "\n";
		}
		$columnSampleCountHash{'Non-contaminant transcript reads'}->{$sample} = $typeCountHash{'non-contaminant'};
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.transcript.remocon.total_read_count.txt", "$directoryPrefix.transcriptome.remocon.total_read_count.txt")) {
		chomp(my $count = `cat $file`);
		print STDERR join("\t", $sample, 'different read count', 'transcript'), "\n" if($count != $columnSampleCountHash{'Non-contaminant transcriptome reads'}->{$sample});
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.transcript.total_read_count.txt", "$directoryPrefix.transcriptome.total_read_count.txt")) {
		chomp(my $count = `cat $file`);
		print STDERR join("\t", $sample, 'different read count', 'transcript'), "\n" if($count != $columnSampleCountHash{'QC-passed reads'}->{$sample});
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.transcript.mapping_strand.read_count.txt", "$directoryPrefix.transcript.remocon.mapping_strand.read_count.txt", "$directoryPrefix.transcriptome.mapping_strand.read_count.txt", "$directoryPrefix.transcriptome.remocon.mapping_strand.read_count.txt")) {
		addColumn('Transcript-mapped reads', '% Transcript-mapped reads', 'Reverse-stranded reads', '% Reverse-stranded reads');
		chomp($columnSampleCountHash{'Transcript-mapped reads'}->{$sample} = `cut -f2 $file | awk '{a += \$o} END {print a}'`);
		$columnSampleCountHash{'% Transcript-mapped reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Transcript-mapped reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
		chomp($columnSampleCountHash{'Reverse-stranded reads'}->{$sample} = `awk -F'\\t' '(\$1 == "reverse") {print \$2}' $file`);
		$columnSampleCountHash{'% Reverse-stranded reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Reverse-stranded reads'}->{$sample} / $columnSampleCountHash{'Transcript-mapped reads'}->{$sample} * 100);
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.transcript.insert_size_metrics.txt", "$directoryPrefix.transcript.remocon.insert_size_metrics.txt", "$directoryPrefix.transcriptome.insert_size_metrics.txt", "$directoryPrefix.transcriptome.remocon.insert_size_metrics.txt")) {
		addColumn('Median insert size', 'Median absolute deviation', 'Mean insert size', 'Standard deviation');
		open(my $reader, "grep -v '^#' $file | grep -v '^\$' |");
		chomp(my $line = <$reader>);
		my @columnList = split(/\t/, $line, -1);
		my %tokenHash = ();
		chomp($line = <$reader>);
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		close($reader);
		$columnSampleCountHash{'Median insert size'}->{$sample} = $tokenHash{'MEDIAN_INSERT_SIZE'};
		$columnSampleCountHash{'Median absolute deviation'}->{$sample} = $tokenHash{'MEDIAN_ABSOLUTE_DEVIATION'};
		$columnSampleCountHash{'Mean insert size'}->{$sample} = sprintf('%.2f', $tokenHash{'MEAN_INSERT_SIZE'});
		$columnSampleCountHash{'Standard deviation'}->{$sample} = sprintf('%.2f', $tokenHash{'STANDARD_DEVIATION'});
	}
	if(-s "$directoryPrefix.rRNA.read_count.txt") {
		addColumn('rRNA-mapped reads', '% rRNA-mapped reads');
		chomp($columnSampleCountHash{'rRNA-mapped reads'}->{$sample} = `cat $directoryPrefix.rRNA.read_count.txt`);
		$columnSampleCountHash{'% rRNA-mapped reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'rRNA-mapped reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.STAR/Log.final.out", "$directoryPrefix.STAR.hg38/Log.final.out")) {
		chomp(my @typeCountList = `cat $file`);
		my %typeCountHash = map {$_->[0] => $_->[1]} grep {defined($_->[1]) && $_->[1] =~ /^[0-9]+$/} map {[split(/\t/, $_, 2)]} @typeCountList;
		my $count = $typeCountHash{'                          Number of input reads |'};
		if($count * 2 == $columnSampleCountHash{'QC-passed reads'}->{$sample}) {
			$typeCountHash{$_} = $typeCountHash{$_} * 2 foreach(keys %typeCountHash);
		} elsif($count != $columnSampleCountHash{'QC-passed reads'}->{$sample}) {
			print STDERR join("\t", $sample, 'different read count', 'genome'), "\n";
		}
		addColumn('Genome-mapped reads', '% Genome-mapped reads', 'Uniquely mapped reads', '% Uniquely mapped reads');
		$columnSampleCountHash{'Genome-mapped reads'}->{$sample} = sum(@typeCountHash{'                   Uniquely mapped reads number |', '        Number of reads mapped to multiple loci |', '        Number of reads mapped to too many loci |'});
		$columnSampleCountHash{'% Genome-mapped reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Genome-mapped reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
		$columnSampleCountHash{'Uniquely mapped reads'}->{$sample} = $typeCountHash{'                   Uniquely mapped reads number |'};
		$columnSampleCountHash{'% Uniquely mapped reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Uniquely mapped reads'}->{$sample} / $columnSampleCountHash{'Genome-mapped reads'}->{$sample} * 100);
	}
	if(-s "$directoryPrefix.genome.total_read_count.txt") {
		chomp(my $count = `cat $directoryPrefix.genome.total_read_count.txt`);
		print STDERR join("\t", $sample, 'different read count', 'genome'), "\n" if($count != $columnSampleCountHash{'QC-passed reads'}->{$sample});
	}
	if(-s "$directoryPrefix.genome.mapped_read_count.txt") {
		addColumn('Genome-mapped reads', '% Genome-mapped reads');
		chomp($columnSampleCountHash{'Genome-mapped reads'}->{$sample} = `cat $directoryPrefix.genome.mapped_read_count.txt`);
		$columnSampleCountHash{'% Genome-mapped reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Genome-mapped reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
	}
	if(-s "$directoryPrefix.genome.uniquely_mapped_read_count.txt") {
		addColumn('Uniquely mapped reads', '% Uniquely mapped reads');
		chomp($columnSampleCountHash{'Uniquely mapped reads'}->{$sample} = `cat $directoryPrefix.genome.uniquely_mapped_read_count.txt`);
		$columnSampleCountHash{'% Uniquely mapped reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Uniquely mapped reads'}->{$sample} / $columnSampleCountHash{'Genome-mapped reads'}->{$sample} * 100);
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.remocon.log", "$directoryPrefix.STAR.remocon.log")) {
		chomp(my @typeCountList = `cat $file`);
		my %typeCountHash = map {$_->[0] => $_->[1]} map {[split(/: /, $_, 2)]} @typeCountList;
		my $count = sum(values %typeCountHash);
		if($count * 2 == $columnSampleCountHash{'QC-passed reads'}->{$sample}) {
			$typeCountHash{$_} = $typeCountHash{$_} * 2 foreach(keys %typeCountHash);
		} elsif($count != $columnSampleCountHash{'QC-passed reads'}->{$sample}) {
			print STDERR join("\t", $sample, 'different read count', 'genome'), "\n";
		}
		addColumn('Contaminant reads', '% Contaminant reads', 'Ambiguous reads', '% Ambiguous reads');
		$columnSampleCountHash{'Contaminant reads'}->{$sample} = $typeCountHash{'contaminant'};
		$columnSampleCountHash{'% Contaminant reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Contaminant reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
		$columnSampleCountHash{'Ambiguous reads'}->{$sample} = $typeCountHash{'ambiguous'};
		$columnSampleCountHash{'% Ambiguous reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Ambiguous reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
	}
	if(my ($file) = grep {-s $_} ("$directoryPrefix.dupMarked.metrics", "$directoryPrefix.remocon.dupMarked.metrics")) {
		my $count = 0;
		my @columnList = ();
		open(my $reader, $file);
		while(my $line = <$reader>) {
			chomp($line);
			if($line eq '') {
				@columnList = ();
			} elsif($line =~ /^## METRICS CLASS\t/) {
				chomp(my $line = <$reader>);
				@columnList = split(/\t/, $line, -1);
			} elsif(@columnList) {
				my %tokenHash = ();
				@tokenHash{@columnList} = split(/\t/, $line, -1);
				$count += $tokenHash{'UNPAIRED_READS_EXAMINED'} + $tokenHash{'READ_PAIRS_EXAMINED'} * 2 + $tokenHash{'UNMAPPED_READS'};
				$columnSampleCountHash{'Mapped reads'}->{$sample} += $tokenHash{'UNPAIRED_READS_EXAMINED'} + $tokenHash{'READ_PAIRS_EXAMINED'} * 2;
				$columnSampleCountHash{'Non-duplicate reads'}->{$sample} += ($tokenHash{'UNPAIRED_READS_EXAMINED'} - $tokenHash{'UNPAIRED_READ_DUPLICATES'}) + ($tokenHash{'READ_PAIRS_EXAMINED'} - $tokenHash{'READ_PAIR_DUPLICATES'}) * 2;
			}
		}
		close($reader);
		if(defined($columnSampleCountHash{'Contaminant reads'}->{$sample})) {
			if($count != $columnSampleCountHash{'QC-passed reads'}->{$sample} - $columnSampleCountHash{'Contaminant reads'}->{$sample}) {
				print STDERR join("\t", $sample, 'different read count', 'genome'), "\n";
			}
			
		} else {
			if($count != $columnSampleCountHash{'QC-passed reads'}->{$sample}) {
				print STDERR join("\t", $sample, 'different read count', 'genome'), "\n";
			}
		}
		addColumn('Mapped reads', '% Mapped reads', 'Non-duplicate reads', '% Non-duplicate reads');
		$columnSampleCountHash{'% Mapped reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Mapped reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
		$columnSampleCountHash{'% Non-duplicate reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Non-duplicate reads'}->{$sample} / $columnSampleCountHash{'QC-passed reads'}->{$sample} * 100);
	}
	if(-s "$directoryPrefix.target.read_count.txt") {
		addColumn('Target reads');
		chomp($columnSampleCountHash{'Target reads'}->{$sample} = `cat $directoryPrefix.target.read_count.txt`);
		if(defined($columnSampleCountHash{'Non-duplicate reads'}->{$sample})) {
			addColumn('% Target reads');
			$columnSampleCountHash{'% Target reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Target reads'}->{$sample} / $columnSampleCountHash{'Non-duplicate reads'}->{$sample} * 100);
		} elsif(defined($columnSampleCountHash{'Mapped reads'}->{$sample})) {
			addColumn('% Target reads');
			$columnSampleCountHash{'% Target reads'}->{$sample} = sprintf('%.1f%%', $columnSampleCountHash{'Target reads'}->{$sample} / $columnSampleCountHash{'Mapped reads'}->{$sample} * 100);
		}
	}
	if(-s "$directoryPrefix.target.depth.txt") {
		chomp(my @typeCountList = `cat $directoryPrefix.target.depth.txt`);
		foreach(map {[split(/\t/, $_, 2)]} @typeCountList) {
			addColumn($_->[0]);
			chomp($columnSampleCountHash{$_->[0]}->{$sample} = $_->[1]);
		}
	}
}

print join("\t", @columnList), "\n";
foreach my $sample (@sampleList) {
	print join("\t", map {defined($_) ? $_ : ''} map {$_->{$sample}} @columnSampleCountHash{@columnList}), "\n";
}

sub addColumn {
	foreach my $column (@_) {
		unless($columnHash{$column}) {
			push(@columnList, $column);
			$columnHash{$column} = 1;
		}
	}
}
