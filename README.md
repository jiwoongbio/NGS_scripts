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

# List gene / gene ID pairs
time perl gff_extract.pl -E */*_genomic.gff.gz gene Dbxref | table.delimitLines.pl - 1 | sed -n 's/GeneID://p' | sort -u > gene.gene_id.txt

# List gene / representative transcript pairs
time perl gff_extract.pl -E */*_genomic.gff.gz transcript_id tag | awk -F'\t' '($2 == "MANE Select")' | table.search.pl - 0 gene.transcript.txt 1 > gene.transcript.MANE_Select.txt
time perl gff_extract.pl -E */*_genomic.gff.gz transcript_id tag | awk -F'\t' '($2 == "RefSeq Select")' | table.search.pl - 0 gene.transcript.txt 1 > gene.transcript.RefSeq_Select.txt
time table.search.pl -v gene.transcript.MANE_Select.txt 0,1 gene.transcript.RefSeq_Select.txt 0,1 | cat gene.transcript.MANE_Select.txt - > gene.transcript.select.txt
```
