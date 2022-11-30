## Requirements

1. Annomen - https://github.com/jiwoongbio/Annomen
2. BWA - http://bio-bwa.sourceforge.net
3. Samtools - http://www.htslib.org
4. STAR - https://github.com/alexdobin/STAR


## Setting for human RNA-seq data analysis (generate_scripts.RNAseq.hg38.pl)

```
git clone https://github.com/jiwoongbio/Annomen.git
mv Annomen Annomen.hg38
cd Annomen.hg38

# Generate genome annotation table
time ./Annomen_table.hg38.sh

# Generate genome index
time bwa index hg38.fasta
time samtools faidx hg38.fasta

# Generate transcriptome index
time bwa index human.rna.fna
time samtools faidx human.rna.fna

# List rRNA transcripts
time sed -n 's/^>//p' human.rna.fna | sed 's/\s/\t/' | sed -r 's/, ([^,]*)$/\t\1/' | awk -F'\t' '($3 == "ribosomal RNA" || $3 == "rRNA")' > rRNA.txt

# List gene / transcript pairs
time table.rearrangeColumns.pl -c Annomen_table.hg38.txt geneName transcriptId | awk '(NR > 1)' | sed -r 's/_([0-9]+\.[0-9]+)/\t\1/' | sort -t $'\t' -k1,1 -k2,2 -k3,3n | uniq | awk -F'\t' -vOFS='\t' '{print $1, $2"_"$3}' > gene.transcript.txt

# Extract transcript sequences with gene annotation
time sed 's/ .*$//' human.rna.fna | tr '\n' '\t' | sed 's/\t$/\n/' | sed 's/\t>/\n>/g' | sed 's/^>//' | table.search.pl gene.transcript.txt 1 - 0 | sed 's/^/>/' | sed 's/\t/\n/g' > transcript.fasta

# Generate transcriptome index
time bwa index transcript.fasta
time samtools faidx transcript.fasta

# Fix GTF file
time gzip -dc GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz | grep -v '^#' | table.substitute_value.pl -i 0 -f chromosome.hg38.txt -o - > genomic.gtf

# Generate STAR index
rm -rf STAR; mkdir STAR; time STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles hg38.fasta --sjdbGTFfile genomic.gtf
```
