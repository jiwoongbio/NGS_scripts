# Usage: Rscript aldex_wilcox.R read_count.tsv output_prefix cutoff_FDR cutoff_effect plot_width plot_height comma_separated_control_samples comma_separated_test_samples

# https://www.bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html

library(ALDEx2)

args <- commandArgs(TRUE)

reads.file <- args[1]
prefix <- args[2]
cutoff.pval <- as.numeric(args[3])
cutoff.effect <- as.numeric(args[4])
width <- as.numeric(args[5])
height <- as.numeric(args[6])
samples.list <- lapply(args[-(1:6)], function(x) {if(file.exists(x)) {scan(x, what = "character", sep = ",")} else {unlist(strsplit(x, ","))}})

reads <- read.table(reads.file, header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(reads) <- reads[, 1]
reads <- reads[, -1]

colnames <- colnames(reads)
colnames(reads) <- 1:length(colnames)

colnumbers <- c()
conditions <- c()
for(index in 1:length(samples.list)) {
	for(colnumber in colnames(reads)[colnames %in% samples.list[[index]]]) {
		colnumbers <- c(colnumbers, colnumber)
		conditions <- c(conditions, index)
	}
}
reads <- reads[, colnumbers]

x <- aldex(reads, conditions, test = "t", effect = TRUE)
x$significant <- ifelse(x$wi.eBH <= cutoff.pval & x$effect >= cutoff.effect, 1, 0)

write.table(cbind(feature = rownames(x), x), file = paste(prefix, "aldex_wilcox.tsv", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

pdf(file = paste(prefix, "aldex_wilcox.MA.pdf", sep = "."), width = width, height = height)
aldex.plot(x, type = "MA", test = "wilcox", cutoff.pval = cutoff.pval, cutoff.effect = cutoff.effect, xlab = "Log-ratio abundance", ylab = "Difference")
dev.off()

pdf(file = paste(prefix, "aldex_wilcox.MW.pdf", sep = "."), width = width, height = height)
aldex.plot(x, type = "MW", test = "wilcox", cutoff.pval = cutoff.pval, cutoff.effect = cutoff.effect, xlab = "Dispersion", ylab = "Difference")
dev.off()
