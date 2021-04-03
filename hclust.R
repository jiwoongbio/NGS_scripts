args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]

table <- data.matrix(table)
if(identical(rownames(table), colnames(table)) && isSymmetric(table)) {
	hc <- hclust(as.dist(table))
} else {
	hc <- hclust(dist(t(table)))
}

pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
plot(hc)
dev.off()
