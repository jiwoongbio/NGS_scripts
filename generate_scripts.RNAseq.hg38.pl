#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
	'p=i' => \(my $threads = 8),
	'P' => \(my $filterPair = ''),
	'T' => \(my $noTranscriptomeMapping = ''),
	'f=f' => \(my $outFilterMinOverLread = ''),
	'n' => \(my $namesort = ''),
	's=s' => \(my $stranded = 'no'),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] sample [...]

Options: -h       display this help message
         -p INT   number of threads [$threads]
         -P       filter pair
         -T       no transcriptome mapping
         -f FLOAT STAR outFilterScoreMinOverLread, outFilterMatchNminOverLread
         -n       HTSeq-count input namesorted BAM file
         -s STR   HTSeq-count stranded [$stranded]

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
#!/usr/bin/bash

EOF
	print $writer <<EOF;
export PATH=/archive/PCDC/shared/jkim23/projects/NGS_scripts:\$PATH

EOF
	print $writer <<EOF;
ANNOMEN=/archive/PCDC/PCDC_Core/shared/pipelines/projects/Annomen.hg38
PICARD=/archive/PCDC/PCDC_Core/shared/pipelines/src/picard/2.27.5/picard.jar

module load trimgalore/0.6.4
module load bwa/intel/0.7.17
module load samtools/intel/1.10
module load star/2.7.3a
module load HTSeq/0.6.1

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

trim_galore --illumina --cores $threads --paired $prefix.pair_filtered.1.fastq.gz $prefix.pair_filtered.2.fastq.gz

mv $prefix.pair_filtered.1_val_1.fq.gz $prefix.1_val_1.fq.gz && rm $prefix.pair_filtered.1.fastq.gz
mv $prefix.pair_filtered.2_val_2.fq.gz $prefix.2_val_2.fq.gz && rm $prefix.pair_filtered.2.fastq.gz

time fastq.read_length.read_count.pl $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz > $prefix.trimgalore.read_length.read_count.txt

EOF
		} else {
			print $writer <<EOF;
trim_galore --illumina --cores $threads --paired $prefix.1.fastq.gz $prefix.2.fastq.gz

time fastq.read_length.read_count.pl $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz > $prefix.trimgalore.read_length.read_count.txt

EOF
		}
		if($noTranscriptomeMapping eq '') {
			print $writer <<EOF;
# Transcriptome mapping for QC
time bwa mem -t $threads \$ANNOMEN/transcript.fasta $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz | samtools sort - > $prefix.transcriptome.sorted.bam && samtools index $prefix.transcriptome.sorted.bam

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
# Genome mapping using STAR
rm -rf $prefix.STAR; mkdir $prefix.STAR; time STAR --runMode alignReads --runThreadN $threads --genomeDir \$ANNOMEN/STAR --readFilesIn $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $prefix.STAR/ --twopassMode Basic --outSAMattributes All --outSAMstrandField intronMotif --outSAMunmapped Within KeepPairs

EOF
		} else {
			print $writer <<EOF;
# Genome mapping using STAR
rm -rf $prefix.STAR; mkdir $prefix.STAR; time STAR --runMode alignReads --runThreadN $threads --genomeDir \$ANNOMEN/STAR --readFilesIn $prefix.1_val_1.fq.gz $prefix.2_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $prefix.STAR/ --twopassMode Basic --outSAMattributes All --outSAMstrandField intronMotif --outSAMunmapped Within KeepPairs --outFilterScoreMinOverLread $outFilterMinOverLread --outFilterMatchNminOverLread $outFilterMinOverLread

EOF
		}
		print $writer <<EOF;
# samtools sort/index
time samtools sort --threads $threads $prefix.STAR/Aligned.out.sam > $prefix.STAR.sorted.bam && samtools index $prefix.STAR.sorted.bam && rm $prefix.STAR/Aligned.out.sam

EOF
	}
	unless(scalar(@prefixList) == 1 && $prefixList[0] eq $sample) {
		print $writer <<EOF;
# Merging BAM files
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
# HTSeq-count
time htseq-count --format=bam --order=pos --stranded=$stranded --type=exon --idattr=gene_id --mode=union $sample.STAR.sorted.bam \$ANNOMEN/genomic.gtf > $sample.htseq-count.txt

EOF
	} else {
		print $writer <<EOF;
time samtools sort -n --threads $threads $sample.STAR.sorted.bam > $sample.STAR.namesorted.bam

# HTSeq-count
time htseq-count --format=bam --order=name --stranded=$stranded --type=exon --idattr=gene_id --mode=union $sample.STAR.namesorted.bam \$ANNOMEN/genomic.gtf > $sample.htseq-count.txt

EOF
	}
	print $writer <<EOF;
time CPM.pl $sample.htseq-count.txt > $sample.CPM.txt
time RPKM.pl \$ANNOMEN/genomic.gtf $sample.htseq-count.txt > $sample.RPKM.txt

EOF
	print $writer <<EOF;
unset ANNOMEN
unset PICARD

module unload trimgalore/0.6.4
module unload bwa/intel/0.7.17
module unload samtools/intel/1.10
module unload star/2.7.3a
module unload HTSeq/0.6.1

EOF
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
