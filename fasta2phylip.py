import sys
(inputFile, outputFile) = (sys.argv[1], sys.argv[2])

from Bio import SeqIO
SeqIO.write(SeqIO.parse(inputFile, "fasta"), outputFile, "phylip")
