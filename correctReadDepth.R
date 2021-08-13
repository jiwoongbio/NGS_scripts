args <- commandArgs(TRUE)

library(TitanCNA)
library(GenomeInfoDb)
library(GenomicRanges)
library(HMMcopy)

source("/archive/PCDC/shared/jkim23/src/TitanCNA/R/utils.R")
#    names(targetIR) <- setGenomeStyle(names(targetIR), genomeStyle)
#    seqnames(targetIR) <- setGenomeStyle(seqnames(targetIR), genomeStyle)

if(length(args) == 6) {
	targetedSequence <- read.delim(args[6], header = FALSE, as.is = TRUE)
	write.table(correctReadDepth(args[1], args[2], args[3], args[4], genomeStyle = "UCSC", targetedSequence = targetedSequence), file = args[5], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
	write.table(correctReadDepth(args[1], args[2], args[3], args[4], genomeStyle = "UCSC"), file = args[5], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
