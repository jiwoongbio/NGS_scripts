#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

my $annomenDirectory = '/archive/PCDC/PCDC_Core/shared/pipelines/projects/Annomen.hg38';
my $picardFile = '/archive/PCDC/PCDC_Core/shared/pipelines/src/picard/2.27.5/picard.jar';
my $gatkFile = '/archive/PCDC/PCDC_Core/shared/pipelines/src/gatk-4.3.0.0/gatk';

my @moduleList = (
	'trimgalore/0.6.4',
	'bwa/intel/0.7.17',
	'samtools/intel/1.10',
	'star/2.7.3a',
	'HTSeq/0.6.1',
);

my $additionalPaths = join(':', $directory);

GetOptions(
	'h' => \(my $help = ''),
	'p=i' => \(my $threads = 8),
	'a=s' => \$annomenDirectory,
	'P' => \(my $filterPair = ''),
	'T' => \(my $noTranscriptomeMapping = ''),
	'f=f' => \(my $outFilterMinOverLread = ''),
	'n' => \(my $namesort = ''),
	's=s' => \(my $stranded = 'no'),
	'v' => \(my $callVariants = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] sample [...]

Options: -h       display this help message
         -p INT   number of threads [$threads]
         -a DIR   Annomen directory [$annomenDirectory]
         -P       filter pair
         -T       no transcriptome mapping
         -f FLOAT STAR outFilterScoreMinOverLread, outFilterMatchNminOverLread
         -n       HTSeq-count input namesorted BAM file
         -s STR   HTSeq-count stranded [$stranded]
         -v       call variants

EOF
}
my @samplePrefixListList = ();
foreach my $sample (@ARGV) {
	if(-d $sample) {
		if(my @prefixList = getPrefixList($sample)) {
			push(@samplePrefixListList, [$sample, @prefixList]);
		} else {
			print STDERR "$sample does not contain fastq files.\n";
		}
	} elsif(-r $sample) {
		open(my $reader, $sample);
		while(my $line = <$reader>) {
			chomp($line);
			my $sample = $line;
			if(-d $sample) {
				if(my @prefixList = getPrefixList($sample)) {
					push(@samplePrefixListList, [$sample, @prefixList]);
				} else {
					print STDERR "$sample does not contain fastq files.\n";
				}
			} else {
				print STDERR "$sample is not a directory.\n";
			}
		}
		close($reader);
	} else {
		print STDERR "$sample is not a directory.\n";
	}
}
foreach(@samplePrefixListList) {
	my ($sample, @prefixList) = @$_;
	open(my $writer, "> $sample/$sample.sh");
	print $writer <<EOF;
#!/bin/bash

EOF
	print $writer <<EOF;
ANNOMEN=$annomenDirectory
PICARD=$picardFile
GATK=$gatkFile

EOF
	print $writer "module load $_\n" foreach(@moduleList);
	print $writer "\n";
	print $writer <<EOF;
export PATH=$additionalPaths:\$PATH

EOF
	print $writer <<EOF;
mkdir -p tmp
export TMPDIR=`readlink -f tmp`

EOF
	foreach my $prefix (@prefixList) {
		print $writer <<EOF;
time fastq.read_length.read_count.pl $prefix.1.fastq.gz $prefix.2.fastq.gz > $prefix.read_length.read_count.txt

EOF
		if($filterPair) {
			print $writer <<EOF;
time fastq.filter_pair.pl $prefix.1.fastq.gz $prefix.2.fastq.gz $prefix.pair_filtered.1.fastq.gz $prefix.pair_filtered.2.fastq.gz > $prefix.pair_filtered.read_pair_count.txt

time fastq.read_length.read_count.pl $prefix.pair_filtered.1.fastq.gz $prefix.pair_filtered.2.fastq.gz > $prefix.pair_filtered.read_length.read_count.txt

# Trim Galore: quality and adapter trimming
trim_galore --illumina --cores $threads --paired $prefix.pair_filtered.1.fastq.gz $prefix.pair_filtered.2.fastq.gz

mv $prefix.pair_filtered.1_val_1.fq.gz $prefix.1_val_1.fq.gz && rm $prefix.pair_filtered.1.fastq.gz
mv $prefix.pair_filtered.2_val_2.fq.gz $prefix.2_val_2.fq.gz && rm $prefix.pair_filtered.2.fastq.gz

time fastq.read_length.read_count.pl $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz > $prefix.trimgalore.read_length.read_count.txt

EOF
		} else {
			print $writer <<EOF;
# Trim Galore: quality and adapter trimming
trim_galore --illumina --cores $threads --paired $prefix.1.fastq.gz $prefix.2.fastq.gz

time fastq.read_length.read_count.pl $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz > $prefix.trimgalore.read_length.read_count.txt

EOF
		}
		if($noTranscriptomeMapping eq '') {
			print $writer <<EOF;
# BWA: mapping to transcriptome
time bwa mem -t $threads \$ANNOMEN/transcript.fasta $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz | samtools sort - > $prefix.transcriptome.sorted.bam && samtools index $prefix.transcriptome.sorted.bam

# Picard: estimate insert size
time java -Djava.io.tmpdir=\$TMPDIR -Xmx16g -jar \$PICARD CollectInsertSizeMetrics INPUT=$prefix.transcriptome.sorted.bam OUTPUT=$prefix.transcriptome.insert_size_metrics.txt HISTOGRAM_FILE=$prefix.transcriptome.insert_size_histogram.pdf VERBOSITY=ERROR QUIET=true VALIDATION_STRINGENCY=LENIENT

time samtools view -c -F 2304 $prefix.transcriptome.sorted.bam > $prefix.transcriptome.total_read_count.txt
time sam.mapping_strand.read_count.pl $prefix.transcriptome.sorted.bam > $prefix.transcriptome.mapping_strand.read_count.txt

time samtools view -f 3 -F 2316 $prefix.transcriptome.sorted.bam | cut -f3 | uniq -c | sed -r 's/^ *([0-9]+) (.*)\$/\\2\\t\\1/' > $prefix.transcriptome.transcript_read_count.txt
time table.search.pl \$ANNOMEN/rRNA.txt 0 $prefix.transcriptome.transcript_read_count.txt 0 | awk -F'\\t' '{a += \$2} END {print a}' > $prefix.rRNA.read_count.txt

rm $prefix.transcriptome.sorted.bam $prefix.transcriptome.sorted.bam.bai

EOF
		}
		if($outFilterMinOverLread eq '') {
			print $writer <<EOF;
# STAR: mapping to genome
rm -rf $prefix.STAR; mkdir $prefix.STAR; time STAR --runMode alignReads --runThreadN $threads --genomeDir \$ANNOMEN/STAR --readFilesIn $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $prefix.STAR/ --twopassMode Basic --outSAMattributes All --outSAMstrandField intronMotif --outSAMunmapped Within KeepPairs

EOF
		} else {
			print $writer <<EOF;
# STAR: mapping to genome
rm -rf $prefix.STAR; mkdir $prefix.STAR; time STAR --runMode alignReads --runThreadN $threads --genomeDir \$ANNOMEN/STAR --readFilesIn $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $prefix.STAR/ --twopassMode Basic --outSAMattributes All --outSAMstrandField intronMotif --outSAMunmapped Within KeepPairs --outFilterScoreMinOverLread $outFilterMinOverLread --outFilterMatchNminOverLread $outFilterMinOverLread

EOF
		}
		print $writer <<EOF;
# Samtools: sort alignment file
time samtools sort --threads $threads $prefix.STAR/Aligned.out.sam > $prefix.STAR.sorted.bam && samtools index $prefix.STAR.sorted.bam && rm $prefix.STAR/Aligned.out.sam

EOF
	}
	unless(scalar(@prefixList) == 1 && $prefixList[0] eq $sample) {
		print $writer <<EOF;
# Samtools: merge BAM files
time samtools merge -f $sample.STAR.sorted.bam `for prefix in @prefixList; do ls \$prefix.STAR.sorted.bam; done` && samtools index $sample.STAR.sorted.bam && for prefix in @prefixList; do rm \$prefix.STAR.sorted.{bam,bam.bai}; done

EOF
		print $writer <<EOF;
for prefix in @prefixList; do cat \$prefix.read_length.read_count.txt; done > $sample.read_length.read_count.txt
for prefix in @prefixList; do cat \$prefix.trimgalore.read_length.read_count.txt; done > $sample.trimgalore.read_length.read_count.txt

for prefix in @prefixList; do awk -vOFS='\\t' '{print NR, \$o}' \$prefix.transcriptome.mapping_strand.read_count.txt; done | table.mergeLines.pl -f sum - 0,1 | cut -f2,3 > $sample.transcriptome.mapping_strand.read_count.txt

for prefix in @prefixList; do cat \$prefix.rRNA.read_count.txt; done | awk '{a += \$o} END {print a}' > $sample.rRNA.read_count.txt

EOF
		print $writer <<EOF;
mkdir $sample.STAR
for prefix in @prefixList; do cat \$prefix.STAR/Log.final.out; done | awk -F'\\t' '(NF == 2 && \$2 ~ /^[0-9]*\$/)' | table.mergeLines.pl -f sum - 0 > $sample.STAR/Log.final.out

EOF
	}
	if($namesort eq '') {
		print $writer <<EOF;
# HTSeq-count: count reads per gene
time htseq-count --format=bam --order=pos --stranded=$stranded --type=exon --idattr=gene_id --mode=union $sample.STAR.sorted.bam \$ANNOMEN/genome.gtf > $sample.htseq-count.txt

EOF
	} else {
		print $writer <<EOF;
# Samtools: sort alignment file by read name
time samtools sort -n --threads $threads $sample.STAR.sorted.bam > $sample.STAR.namesorted.bam

# HTSeq-count: count reads per gene
time htseq-count --format=bam --order=name --stranded=$stranded --type=exon --idattr=gene_id --mode=union $sample.STAR.namesorted.bam \$ANNOMEN/genome.gtf > $sample.htseq-count.txt

EOF
	}
	print $writer <<EOF;
time CPM.pl $sample.htseq-count.txt > $sample.CPM.txt
time RPKM.pl \$ANNOMEN/genome.gtf $sample.htseq-count.txt > $sample.RPKM.txt

EOF
	if($callVariants) {
		print $writer <<EOF;
time samtools view -b -F 2304 $sample.STAR.sorted.bam > $sample.STAR.sorted.filtered.bam

# Picard: add read group information
time java -Djava.io.tmpdir=\$TMPDIR -Xmx16g -jar \$PICARD AddOrReplaceReadGroups INPUT=$sample.STAR.sorted.filtered.bam OUTPUT=$sample.rgAdded.bam SORT_ORDER=coordinate RGID=$sample RGLB=$sample RGPL=illumina RGPU=$sample RGSM=$sample VERBOSITY=ERROR QUIET=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true && rm $sample.STAR.sorted.filtered.bam

EOF
	print $writer <<EOF;
# Picard: mark PCR duplicates
time java -Djava.io.tmpdir=\$TMPDIR -Xmx16g -jar \$PICARD MarkDuplicates INPUT=$sample.rgAdded.bam OUTPUT=$sample.dupMarked.bam METRICS_FILE=$sample.dupMarked.metrics VERBOSITY=ERROR QUIET=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true && rm $sample.rgAdded.{bam,bai}

EOF
	print $writer <<EOF;
# GATK: split reads with N in cigar
time \$GATK --java-options "-Djava.io.tmpdir=\$TMPDIR -Xmx16g" SplitNCigarReads -R \$ANNOMEN/genome.fasta -I $sample.dupMarked.bam -O $sample.SplitNCigarReads.bam --verbosity ERROR && rm $sample.dupMarked.{bam,bai}

EOF
	print $writer <<EOF;
# GATK: base quality score recalibration (BQSR)
time \$GATK --java-options "-Djava.io.tmpdir=\$TMPDIR -Xmx16g" BaseRecalibrator -R \$ANNOMEN/genome.fasta -I $sample.SplitNCigarReads.bam -O $sample.recal_data.table --use-original-qualities --known-sites \$ANNOMEN/snp_b155.vcf.gz --verbosity ERROR
time \$GATK --java-options "-Djava.io.tmpdir=\$TMPDIR -Xmx16g" ApplyBQSR -R \$ANNOMEN/genome.fasta -I $sample.SplitNCigarReads.bam --use-original-qualities --bqsr-recal-file $sample.recal_data.table -O $sample.recal_reads.bam --verbosity ERROR && rm $sample.SplitNCigarReads.{bam,bai}

# GATK: call variants with HaplotypeCaller
time \$GATK --java-options "-Djava.io.tmpdir=\$TMPDIR -Xmx16g" HaplotypeCaller -R \$ANNOMEN/genome.fasta -I $sample.recal_reads.bam -O $sample.raw.snps.indels.g.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
time \$GATK --java-options "-Djava.io.tmpdir=\$TMPDIR -Xmx16g" GenotypeGVCFs -R \$ANNOMEN/genome.fasta -V $sample.raw.snps.indels.g.vcf.gz -O $sample.raw.snps.indels.vcf.gz -D \$ANNOMEN/snp_b155.vcf.gz

EOF
	print $writer <<EOF;
# Custom Perl script: variant filtering with cutoffs (QD < 2, FS > 60, MQ < 40, DP < 3, GQ < 7)
time perl \$ANNOMEN/vcf.filter.pl $sample.raw.snps.indels.vcf.gz QD 'QD < 2' \\
	| perl \$ANNOMEN/vcf.filter.pl - FS 'FS > 60' \\
	| perl \$ANNOMEN/vcf.filter.pl - MQ 'MQ < 40' \\
	| perl \$ANNOMEN/vcf.filter.pl -g - DP 'DP < 3' \\
	| perl \$ANNOMEN/vcf.filter.pl -g - GQ 'GQ < 7' \\
	> $sample.filtered.vcf

EOF
	print $writer <<EOF;
# Custom Perl script: variant annotation (RefSeq genes, tandem repeat)
time perl \$ANNOMEN/leftalignIndel.pl $sample.filtered.vcf \$ANNOMEN/genome.fasta | perl \$ANNOMEN/sort_by_reference.pl - \$ANNOMEN/genome.fasta 0 1 \\
	| perl \$ANNOMEN/Annomen.pl - \$ANNOMEN/genome.fasta \$ANNOMEN/Annomen_table.txt \$ANNOMEN/refseq.transcript.fasta \$ANNOMEN/refseq.protein.fasta \\
	> $sample.annotated.vcf

# Custom Perl script: generate variant tables
time perl \$ANNOMEN/vcf.table.pl -c \$ANNOMEN/column.name.depth.dbSNP.txt $sample.annotated.vcf > $sample.table.variant.txt

time table.filter.pl $sample.table.variant.txt 'dbSNP' 'eq' '' > $sample.table.variant.not_dbSNP.txt
time table.variant_count.pl $sample.table.variant.txt > $sample.table.variant_count.txt
time table.variant_count.pl $sample.table.variant.not_dbSNP.txt > $sample.table.variant_count.not_dbSNP.txt

EOF
	}
	print $writer <<EOF;
unset ANNOMEN
unset PICARD
unset GATK

EOF
	print $writer "module unload $_\n" foreach(@moduleList);
	print $writer "\n";
	close($writer);
}

sub getPrefixList {
	my ($sample) = @_;
	my @prefixList = ();
	if(-r "$sample/prefix.txt") {
		chomp(@prefixList = `cat $sample/prefix.txt`);
	} else {
		push(@prefixList, $sample);
	}
	my @availablePrefixList = ();
	foreach my $prefix (@prefixList) {
		unless(-r (my $fastqFile = "$sample/$prefix.1.fastq.gz")) {
			print STDERR "Can't find '$fastqFile'\n";
			next;
		}
		unless(-r (my $fastqFile = "$sample/$prefix.2.fastq.gz")) {
			print STDERR "Can't find '$fastqFile'\n";
			next;
		}
		if(scalar(@prefixList) > 1 && $prefix eq $sample) {
			print STDERR "Change prefix '$sample/$prefix' to another.\n";
			next;
		}
		push(@availablePrefixList, $prefix);
	}
	return @availablePrefixList;
}
