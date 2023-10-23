args <- commandArgs(TRUE)

fasta_file <- args[1]
threads <- as.numeric(args[2])

tree_file <- args[3]

library(ips)

exec <- "/archive/PCDC/PCDC_Core/shared/pipelines/src/standard-RAxML-8.2.12/raxmlHPC-PTHREADS"

DNAbin <- read.dna(fasta_file, format = "fasta")
tr <- raxml(DNAbin, m = "GTRCAT", f = "d", N = 2, exec = exec, threads = threads)

write.tree(tr$bestTree, file = tree_file)
