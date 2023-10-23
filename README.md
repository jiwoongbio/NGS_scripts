# NGS_scripts

Scripts for analyzing Next-Generation Sequencing data


## Requirements

1. Perl - https://www.perl.org
2. R - http://www.r-project.org
3. Perl modules "Bio::DB::Fasta", "Bio::DB::Taxonomy" - https://bioperl.org
4. wget - https://www.gnu.org/software/wget/
5. Linux commands: sort, gzip, ...


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.

```
git clone https://github.com/jiwoongbio/NGS_scripts.git
```


## Human RNA-seq data analysis

1. Setting reference genome data

```
# Requirements
# 1. Annomen - https://github.com/jiwoongbio/Annomen
# 2. BWA - http://bio-bwa.sourceforge.net
# 3. Samtools - http://www.htslib.org
# 4. STAR - https://github.com/alexdobin/STAR
#
# - You can use conda to install the requirements as follows:
# conda create -n Annomen -c bioconda perl perl-bioperl emboss
# conda install -n Annomen -c anaconda wget
# conda install -n Annomen -c bioconda bwa samtools star

git clone https://github.com/jiwoongbio/Annomen.git
mv Annomen Annomen.hg38
cd Annomen.hg38

# Generate genome annotation table
time ./Annomen_table.hg38.sh

#PICARD=where/is/picard.jar
time java -jar $PICARD CreateSequenceDictionary REFERENCE=genome.fasta

# Generate genome index
time bwa index genome.fasta
time samtools faidx genome.fasta

# Calculate genome sequence lengths
time fasta.length.pl genome.fasta > genome.length.txt

# Calculate refseq transcript sequence lengths
time fasta.length.pl refseq.transcript.fasta > refseq.transcript.length.txt

# List gene / transcript pairs
time perl gff_extract.pl -E */*_genomic.gff.gz gene transcript_id | table.search.pl refseq.transcript.length.txt 0 - 1 | sed -r 's/_([0-9]+\.[0-9]+)/\t\1/' | sort -t $'\t' -k1,1 -k2,2 -k3,3n | uniq | awk -F'\t' -vOFS='\t' '{print $1, $2"_"$3}' > gene.transcript.txt

# Extract transcript sequences with gene annotation
time sed 's/\t/ /g' refseq.transcript.fasta | tr '\n' '\t' | sed 's/\t$/\n/' | sed 's/\t>/\n>/g' | sed 's/^>//' | sed 's/ /\t/' | table.search.pl gene.transcript.txt 1 - 0 | sed 's/\t/ /' | sed 's/^/>/' | sed 's/\t/\n/g' > transcript.fasta

# Generate transcript index
time bwa index transcript.fasta
time samtools faidx transcript.fasta

# Calculate transcript sequence lengths
time fasta.length.pl transcript.fasta > transcript.length.txt

# List rRNA
time sed 's/\t/ /g' transcript.fasta | sed -n 's/^>//p' | sed 's/ /\t/' | sed -r 's/, ([^,]*)$/\t\1/' | awk -F'\t' '($3 == "ribosomal RNA" || $3 == "rRNA")' > rRNA.txt

# Fix GTF file
time gzip -dc */*_genomic.gtf.gz | grep -v '^#' | table.substitute_value.pl -i 0 -f chromosome.UCSC.txt -o - | table.search.pl genome.length.txt 0 - 0 > genome.gtf

# Generate STAR index
rm -rf STAR; mkdir STAR; time STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles genome.fasta --sjdbGTFfile genome.gtf

# SpliceFisher
time git clone https://github.com/jiwoongbio/SpliceFisher.git
rm SpliceFisher/*.bam
rm SpliceFisher/*.txt

time (cd SpliceFisher; ./prepare.sh ../genome.gtf)
time (cd SpliceFisher; perl region_type.pl ../genome.gtf > region_type.txt)
time (cd SpliceFisher; perl region_type.pl -S forward ../genome.gtf | awk -F'\t' -vOFS='\t' '{print $1, $2, $3, "+", $4}' > region_type.forward.txt)
time (cd SpliceFisher; perl region_type.pl -S reverse ../genome.gtf | awk -F'\t' -vOFS='\t' '{print $1, $2, $3, "-", $4}' > region_type.reverse.txt)

# Generate gene + transcript index
time cat genome.fasta transcript.fasta > genome_transcript.fasta
time bwa index genome_transcript.fasta
time samtools faidx genome_transcript.fasta

# List gene / gene ID pairs
time perl gff_extract.pl -E */*_genomic.gff.gz gene Dbxref | table.delimitLines.pl - 1 | sed -n 's/GeneID://p' | sort -u > gene.gene_id.txt

# List gene / representative transcript pairs
time perl gff_extract.pl -E */*_genomic.gff.gz transcript_id tag | awk -F'\t' '($2 == "MANE Select")' | table.search.pl - 0 gene.transcript.txt 1 > gene.transcript.MANE_Select.txt
time perl gff_extract.pl -E */*_genomic.gff.gz transcript_id tag | awk -F'\t' '($2 == "RefSeq Select")' | table.search.pl - 0 gene.transcript.txt 1 > gene.transcript.RefSeq_Select.txt

time perl gff_extract.pl -E */*_genomic.gff.gz transcript_id tag | awk -F'\t' '($2 == "MANE Plus Clinical")' | table.search.pl - 0 gene.transcript.txt 1 > gene.transcript.MANE_Plus_Clinical.txt
time perl gff_extract.pl -E */*_genomic.gff.gz transcript_id tag | awk -F'\t' '($2 == "RefSeq Plus Clinical")' | table.search.pl - 0 gene.transcript.txt 1 > gene.transcript.RefSeq_Plus_Clinical.txt

time table.search.pl -v gene.transcript.MANE_Select.txt 0,1 gene.transcript.RefSeq_Select.txt 0,1 gene.transcript.MANE_Plus_Clinical.txt 0,1 gene.transcript.RefSeq_Plus_Clinical.txt 0,1 | cat gene.transcript.MANE_Select.txt - > gene.transcript.select.txt

# dbSNP 155
time wget --no-verbose --no-check-certificate https://ftp.ncbi.nlm.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz
time gzip -dc GCF_000001405.39.gz | grep -v '^#' | table.substitute_value.pl -i 0 -f chromosome.UCSC.txt -o - | table.search.pl genome.length.txt 0 - 0 | bash -c "cat <(gzip -dc GCF_000001405.39.gz | head -n1000 | grep '^#') -" | perl leftalignIndel.pl - genome.fasta | perl sort_by_reference.pl - genome.fasta 0 1 | bgzip > snp_b155.vcf.gz
time tabix --preset vcf snp_b155.vcf.gz

# dbSNP 156
time wget --no-verbose --no-check-certificate https://ftp.ncbi.nlm.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz
time gzip -dc GCF_000001405.40.gz | grep -v '^#' | table.substitute_value.pl -i 0 -f chromosome.UCSC.txt -o - | table.search.pl genome.length.txt 0 - 0 | bash -c "cat <(gzip -dc GCF_000001405.40.gz | head -n1000 | grep '^#') -" | perl leftalignIndel.pl - genome.fasta | perl sort_by_reference.pl - genome.fasta 0 1 | bgzip > snp_b156.vcf.gz
time tabix --preset vcf snp_b156.vcf.gz
```

2. Analyzing a RNA-seq sample

```
mkdir SAMPLE
ln -sf RAW_DATA/SAMPLE_R1.fastq.gz SAMPLE/SAMPLE.1.fastq.gz
ln -sf RAW_DATA/SAMPLE_R2.fastq.gz SAMPLE/SAMPLE.2.fastq.gz

generate_scripts.RNAseq.pl SAMPLE

sbatch.pl -p PARTITION_OF_32GB_FREE_RAM SAMPLE
```
