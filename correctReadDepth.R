args <- commandArgs(TRUE)

library(TitanCNA)
if(length(args) == 6) {
	exons <- read.delim(args[6], header = FALSE, as.is = TRUE)
	write.table(correctReadDepth(args[1], args[2], args[3], args[4], genomeStyle = "UCSC", targetedSequence = exons), file = args[5], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
	write.table(correctReadDepth(args[1], args[2], args[3], args[4], genomeStyle = "UCSC"), file = args[5], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
